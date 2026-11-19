// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file reduced_model.c
 * @brief Branch-and-cut on the reduced model (see models.h).
 */
#include "exact/models.h"

#include "common/timer.h"
#include "common/util.h"
#include "exact/grb.h"
#include "exact/plan.h"
#include "exact/routes.h"
#include "exact/separation.h"
#include "exact/subproblem_models.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** A column x_rtv is selected when its value exceeds this threshold. */
static const double SELECTED = 0.1;
/** Time limits of the successive solver rounds (s); the last one is repeated. */
static const double ROUND_LIMITS[] = {1500, 1500, 1500, 1500, 1500, 15000, 15000, 20000, 20000, 1e100};

/** @brief Column x_rtv of the reduced model. */
typedef struct Column {
    long route;
    int8_t period;
    int8_t type;
} Column;

/** @brief A cut stored for the following rounds. */
typedef struct StoredCut {
    int len;
    int *ind;
    double *val;
    double rhs;
} StoredCut;

/** @brief Model layout and callback state. */
typedef struct Reduced {
    ExactContext X;
    GRBenv *env;
    Column *col; /* [nx] */
    int nx;      /* number of x columns */
    int y0;      /* column of y_0 */
    int *z0;     /* [num_installations] column of z_i0 */
    int c0;      /* column of c_{0,0,0}; c_vtl = c0 + (v*P + t)*P + l */
    int ncols;
    int *class_first;       /* [T*P*P + 1] first entry in class_col of each class (type, day, span) */
    int *class_col;         /* [nx] the x columns, grouped by class */
    struct WeekModel *week; /* [T] weights of the classes of each type, built when first needed */
    long week_cuts, pattern_cuts;
    double *sol;  /* [ncols] callback buffer */
    Plan plan;    /* plan being checked */
    Plan reduced; /* scratch plan of Reduce */
    Plan best;    /* best accepted plan */
    double best_obj;
    StoredCut *cuts; /* cuts to add before the next round */
    int num_cuts, cap_cuts;
    GrbRow row;
    /* output while solving, so that a run stopped from outside (memory, job limit) leaves its results */
    const char *output; /* solution file, rewritten for every accepted plan (nullptr: none) */
    bool progress;      /* print a progress line at most every PROGRESS_SECONDS */
    double wall0;       /* start of the run */
    double lb;          /* best bound so far, over the root relaxation and the rounds */
    double next_progress;
    int round;
    bool relaxing; /* root_relaxation() is running: the callback does nothing */
} Reduced;

/** Seconds between two progress lines. */
static const double PROGRESS_SECONDS = 60;

/* Column of c_vtl for vessel type @p type, day @p period and span @p span (spans beyond the week are capped). */
static int c_var(const Reduced *R, int type, int period, int span)
{
    const int P = R->X.inst->num_periods;
    return R->c0 + (type * P + period) * P + imin(span, P - 1);
}

/* ------------------------------------------------------------------------- */
/*  Model construction                                                       */
/* ------------------------------------------------------------------------- */

/* Rotating every pattern by one day must give another pattern of the same installation. */
static bool patterns_rotation_closed(const Instance *in)
{
    const int P = in->num_periods;
    for (int i = 0; i < in->num_installations; i++) {
        const Installation *u = &in->inst[i];
        for (int a = 0; a < u->num_patterns; a++) {
            bool found = false;
            for (int b = 0; b < u->num_patterns && !found; b++) {
                bool same = true;
                for (int t = 0; t < P && same; t++)
                    same = u->pattern[b * P + (t + 1) % P] == u->pattern[a * P + t];
                found = same;
            }
            if (!found)
                return false;
        }
    }
    return true;
}

/* The horizon is cyclic, so any solution can be rotated until a weekly visit of a
 * frequency-1 installation falls on day 1. Returns that installation, or -1 when
 * the symmetry cannot be broken this way. The choice (lowest frequency, then
 * lowest demand) follows the original implementation. */
static int installation_fixed_on_day_one(const Instance *in)
{
    int k = 0;
    for (int i = in->num_installations - 1; i > 0; --i) {
        const Installation *a = &in->inst[i], *b = &in->inst[k];
        if (a->frequency < b->frequency || (a->frequency == b->frequency && a->demand < b->demand))
            k = i;
    }
    if (in->inst[k].frequency != 1 || !patterns_rotation_closed(in))
        return -1;
    return k;
}

/* Lists the x_rtv columns: every (voyage, day, type) with a type able to carry the voyage, except
 * those removed by the symmetry fixing of one installation and by the optional usage filter. */
static void build_columns(Reduced *R)
{
    const Instance *in = R->X.inst;
    const RoutePool *pool = R->X.pool;
    const int fixed = installation_fixed_on_day_one(in);
    const bool filter = in->num_clusters == 1 && R->X.par->min_usage > 0;
    long cap = 1024;
    R->col = XMALLOC_ARRAY(Column, cap);
    R->nx = 0;
    for (long r = 0; r < pool->count; r++) {
        const Route *route = &pool->route[r];
        const bool has_fixed = fixed >= 0 && route_visits(pool, r, fixed);
        for (int t = 0; t < in->num_periods; t++)
            for (int k = 0; k < in->num_types; k++) {
                if (in->type[k].max_vessels == 0 || !route_fits(route, k))
                    continue;
                if (has_fixed && t != 0)
                    continue;
                if (filter && 100.0f * (float)route->load / (float)in->type[k].capacity < R->X.par->min_usage)
                    continue;
                if (R->nx == cap)
                    R->col = xrealloc(R->col, sizeof(Column) * (size_t)(cap *= 2));
                if (R->nx == INT32_MAX)
                    die("the reduced model has too many columns");
                R->col[R->nx++] = (Column){.route = r, .period = (int8_t)t, .type = (int8_t)k};
            }
    }
}

/* Builds the reduced model (23)-(31) over the columns of @p R (variables x, y, z, c, constraints
 * (24)-(27) and the max constraints defining c) and records the column layout in @p R. */
static GRBmodel *build_model(Reduced *R)
{
    const Instance *in = R->X.inst;
    const RoutePool *pool = R->X.pool;
    const int P = in->num_periods, T = in->num_types, N = in->num_installations;
    GRBenv *env = R->env;
    GRBmodel *m = grb_new_model(env, "reduced", R->X.par->threads, R->X.par->solver_log);
    GRBenv *menv = GRBgetenv(m);
    GRB_CHECK(env, GRBsetdblparam(menv, GRB_DBL_PAR_HEURISTICS, 0.01));
    GRB_CHECK(env, GRBsetintparam(menv, GRB_INT_PAR_MIPFOCUS, 3));
    GRB_CHECK(env, GRBsetintparam(menv, GRB_INT_PAR_CUTS, 3));
    GRB_CHECK(env, GRBsetintparam(menv, GRB_INT_PAR_LAZYCONSTRAINTS, 1));

    /* columns: x, y, z, c */
    for (int j = 0; j < R->nx; j++) {
        const Column *c = &R->col[j];
        const double obj = (double)in->type[c->type].var_cost * pool->route[c->route].dist;
        GRB_CHECK(env, GRBaddvar(m, 0, nullptr, nullptr, obj, 0.0, 1.0, GRB_BINARY, nullptr));
    }
    R->y0 = R->nx;
    for (int k = 0; k < T; k++)
        GRB_CHECK(env, GRBaddvar(m, 0, nullptr, nullptr, in->type[k].fixed_cost, 0.0, in->type[k].max_vessels,
                                 GRB_INTEGER, nullptr));
    R->z0 = XMALLOC_ARRAY(int, N);
    int next = R->y0 + T;
    for (int i = 0; i < N; i++) {
        R->z0[i] = next;
        for (int p = 0; p < in->inst[i].num_patterns; p++, next++)
            GRB_CHECK(env, GRBaddvar(m, 0, nullptr, nullptr, 0.0, 0.0, 1.0, GRB_BINARY, nullptr));
    }
    R->c0 = next;
    for (int j = 0; j < T * P * P; j++, next++)
        GRB_CHECK(env, GRBaddvar(m, 0, nullptr, nullptr, 0.0, 0.0, 1.0, GRB_BINARY, nullptr));
    R->ncols = next;
    GRB_CHECK(env, GRBupdatemodel(m));

    /* Rows indexed by (installation, day), (type, day) and (type, day, span), filled column-wise. */
    GrbRow *visit = XCALLOC_ARRAY(GrbRow, N * P), *busy = XCALLOC_ARRAY(GrbRow, T * P), *dep = XCALLOC_ARRAY(GrbRow, P),
           *span = XCALLOC_ARRAY(GrbRow, T * P * P);
    for (int j = 0; j < R->nx; j++) {
        const Column *c = &R->col[j];
        const Route *route = &pool->route[c->route];
        const int *s = route_stops(pool, c->route);
        for (int a = 0; a < route->size; a++)
            row_add(&visit[s[a] * P + c->period], j, 1.0);
        for (int d = 0; d < P; d++)
            if (route_active(in, route, c->period, d))
                row_add(&busy[c->type * P + d], j, 1.0);
        row_add(&dep[c->period], j, 1.0);
        row_add(&span[(c->type * P + c->period) * P + imin(route->span, P - 1)], j, 1.0);
    }
    GrbRow row = {0};
    for (int i = 0; i < N; i++) { /* (24) one pattern per installation */
        row_clear(&row);
        for (int p = 0; p < in->inst[i].num_patterns; p++)
            row_add(&row, R->z0[i] + p, 1.0);
        grb_add_row(m, &row, GRB_EQUAL, 1.0, nullptr);
    }
    for (int i = 0; i < N; i++) /* (25) visits follow the pattern */
        for (int t = 0; t < P; t++) {
            GrbRow *v = &visit[i * P + t];
            for (int p = 0; p < in->inst[i].num_patterns; p++)
                if (inst_pattern_visits(&in->inst[i], P, p, t))
                    row_add(v, R->z0[i] + p, -1.0);
            grb_add_row(m, v, GRB_EQUAL, 0.0, nullptr);
        }
    for (int k = 0; k < T; k++) /* (26) fleet >= voyages active on each day */
        for (int d = 0; d < P; d++) {
            GrbRow *v = &busy[k * P + d];
            row_add(v, R->y0 + k, -1.0);
            grb_add_row(m, v, GRB_LESS_EQUAL, 0.0, nullptr);
        }
    for (int t = 0; t < P; t++) /* (27) departures per day <= D = 2 |N^B| */
        grb_add_row(m, &dep[t], GRB_LESS_EQUAL, (double)MAX_DEPARTURES_PER_BERTH * in->num_berths, nullptr);
    for (int k = 0; k < T; k++) /* c_vtl = max of the x of type v, day t and span l */
        for (int t = 0; t < P; t++)
            for (int l = 0; l < P; l++) {
                const GrbRow *v = &span[(k * P + t) * P + l];
                GRB_CHECK(env, GRBaddgenconstrMax(m, nullptr, c_var(R, k, t, l), v->len, v->ind, 0.0));
            }
    /* the columns of every class, for the cuts on the vessels */
    R->class_first = XMALLOC_ARRAY(int, T * P * P + 1);
    R->class_col = XMALLOC_ARRAY(int, imax(R->nx, 1));
    R->class_first[0] = 0;
    for (int c = 0; c < T * P * P; c++) {
        if (span[c].len > 0)
            memcpy(&R->class_col[R->class_first[c]], span[c].ind, sizeof(int) * (size_t)span[c].len);
        R->class_first[c + 1] = R->class_first[c] + span[c].len;
    }
    for (int j = 0; j < N * P; j++)
        row_free(&visit[j]);
    for (int j = 0; j < T * P; j++)
        row_free(&busy[j]);
    for (int j = 0; j < P; j++)
        row_free(&dep[j]);
    for (int j = 0; j < T * P * P; j++)
        row_free(&span[j]);
    free(visit);
    free(busy);
    free(dep);
    free(span);
    row_free(&row);
    GRB_CHECK(env, GRBupdatemodel(m));
    return m;
}

/* ------------------------------------------------------------------------- */
/*  Callback: separation procedure (Algorithm 1)                             */
/* ------------------------------------------------------------------------- */

/* Adds row <= rhs as a lazy constraint now and stores it for the next rounds. */
static void add_cut(Reduced *R, void *cbdata, const GrbRow *row, double rhs)
{
    GRB_CHECK(R->env, GRBcblazy(cbdata, row->len, row->ind, row->val, GRB_LESS_EQUAL, rhs));
    if (R->num_cuts == R->cap_cuts) {
        R->cap_cuts = R->cap_cuts ? 2 * R->cap_cuts : 64;
        R->cuts = xrealloc(R->cuts, sizeof(StoredCut) * (size_t)R->cap_cuts);
    }
    StoredCut *c = &R->cuts[R->num_cuts++];
    c->len = row->len;
    c->rhs = rhs;
    c->ind = XMALLOC_ARRAY(int, imax(row->len, 1));
    c->val = XMALLOC_ARRAY(double, imax(row->len, 1));
    memcpy(c->ind, row->ind, sizeof(int) * (size_t)row->len);
    memcpy(c->val, row->val, sizeof(double) * (size_t)row->len);
}

/* Cut (32): the voyages of the (conflicting part of the) plan cannot all be selected. */
static void add_nogood(Reduced *R, void *cbdata, const Plan *plan, bool conflicts_only)
{
    GrbRow *row = &R->row;
    row_clear(row);
    for (int i = 0; i < plan->count; i++)
        if (!conflicts_only || plan->trip[i].conflict)
            row_add(row, plan->trip[i].var, 1.0);
    if (row->len == 0) /* no conflict marked: forbid the whole plan, which is always valid */
        for (int i = 0; i < plan->count; i++)
            row_add(row, plan->trip[i].var, 1.0);
    add_cut(R, cbdata, row, row->len - 1.0);
}

/* ---- Cuts on the vessels of one type ------------------------------------------------------------
 *
 * Call class (t, l) the voyages of one type that depart on day t and span l days: they occupy a
 * vessel on the same days. In a week a vessel performs voyages of pairwise disjoint classes. Give
 * every class c a weight w_c >= 0 such that no such week weighs more than 1; then
 *
 *     sum_c w_c * (number of voyages of class c)  <=  y_v,
 *
 * every vessel carrying a weight of at most 1. Constraint (26) is the case w = 1 on the classes
 * active on one day. The weights that make the left-hand side largest for the incumbent are given
 * by a small linear programme, with one variable per class and one row per maximal week; its value
 * is the number of vessels that would be needed if a voyage could be shared between vessels (the
 * fractional chromatic number). When it exceeds y_v the inequality cuts the incumbent off;
 * otherwise the corrected cut (33) below is used. */

/* Do the voyages of classes (day a, span la) and (day b, span lb) share a day? */
static bool classes_overlap(int P, int a, int la, int b, int lb)
{
    return (b - a + P) % P <= la || (a - b + P) % P <= lb;
}

/** Above this number of maximal weeks the linear programme is not built (long horizons). */
enum { MAX_WEEKS = 50000 };

/* The classes of one vessel type and the linear programme over their weights. */
typedef struct WeekModel {
    bool built;
    int m;        /* classes with voyages in the model; -1 when the programme is not used */
    int *cls;     /* [m] class: day * P + span */
    GRBmodel *lp; /* maximise the weight of the incumbent; a row per maximal week */
    int rows;
    int *row_first, *row_class; /* the maximal weeks: classes (indices into cls) of each row */
    int cap_first, cap_class;
} WeekModel;

/* Enumerates the sets of pairwise disjoint classes that contain @p chosen and only add classes
 * from index @p from on, and records the maximal ones. Returns false beyond MAX_WEEKS. */
static bool week_search(WeekModel *w, int P, int from, int *chosen, int size)
{
    for (int i = from; i < w->m; i++) {
        bool is_free = true;
        for (int k = 0; k < size && is_free; k++)
            is_free = !classes_overlap(P, w->cls[i] / P, w->cls[i] % P, w->cls[chosen[k]] / P, w->cls[chosen[k]] % P);
        if (!is_free)
            continue;
        chosen[size] = i;
        if (!week_search(w, P, i + 1, chosen, size + 1))
            return false;
    }
    for (int i = 0; i < w->m; i++) { /* maximal: every other class overlaps one of the week */
        bool is_free = true;
        for (int k = 0; k < size && is_free; k++)
            is_free = chosen[k] != i &&
                      !classes_overlap(P, w->cls[i] / P, w->cls[i] % P, w->cls[chosen[k]] / P, w->cls[chosen[k]] % P);
        if (is_free)
            return true;
    }
    if (w->rows == MAX_WEEKS)
        return false;
    if (w->rows + 2 > w->cap_first) {
        w->cap_first = imax(2 * w->cap_first, 256);
        w->row_first = xrealloc(w->row_first, sizeof(int) * (size_t)w->cap_first);
    }
    if (w->rows == 0)
        w->row_first[0] = 0;
    const int at = w->row_first[w->rows];
    if (at + size > w->cap_class) {
        w->cap_class = imax(2 * w->cap_class, 1024);
        w->row_class = xrealloc(w->row_class, sizeof(int) * (size_t)w->cap_class);
    }
    memcpy(&w->row_class[at], chosen, sizeof(int) * (size_t)size);
    w->row_first[++w->rows] = at + size;
    return true;
}

/* Lists the classes of @p type and their maximal weeks, and builds the linear programme. */
static void week_model_build(Reduced *R, int type)
{
    const int P = R->X.inst->num_periods;
    WeekModel *w = &R->week[type];
    w->built = true;
    w->cls = XMALLOC_ARRAY(int, P * P);
    for (int c = 0; c < P * P; c++)
        if (R->class_first[type * P * P + c] < R->class_first[type * P * P + c + 1])
            w->cls[w->m++] = c;
    int *chosen = XMALLOC_ARRAY(int, imax(w->m, 1));
    const bool complete = week_search(w, P, 0, chosen, 0);
    free(chosen);
    if (!complete || w->m == 0) {
        w->m = -1;
        return;
    }
    w->lp = grb_new_model(R->env, "weeks", 1, false);
    GRB_CHECK(R->env, GRBsetintattr(w->lp, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE));
    for (int i = 0; i < w->m; i++)
        GRB_CHECK(R->env, GRBaddvar(w->lp, 0, nullptr, nullptr, 0.0, 0.0, 1.0, GRB_CONTINUOUS, nullptr));
    GRB_CHECK(R->env, GRBupdatemodel(w->lp));
    GrbRow row = {0};
    for (int r = 0; r < w->rows; r++) {
        row_clear(&row);
        for (int i = w->row_first[r]; i < w->row_first[r + 1]; i++)
            row_add(&row, w->row_class[i], 1.0);
        grb_add_row(w->lp, &row, GRB_LESS_EQUAL, 1.0, nullptr);
    }
    row_free(&row);
}

/* The weighted cut for the voyages of @p type of the plan, if they need more than the fleet even
 * when shared between vessels. Returns true if the cut was added. */
static bool add_week_cut(Reduced *R, void *cbdata, const Plan *plan, int type)
{
    const int P = R->X.inst->num_periods;
    WeekModel *w = &R->week[type];
    if (!w->built)
        week_model_build(R, type);
    if (w->m < 0)
        return false;
    int *count = XCALLOC_ARRAY(int, P * P);
    double *weight = XMALLOC_ARRAY(double, w->m);
    for (int i = 0; i < plan->count; i++) {
        const Trip *t = &plan->trip[i];
        if (t->type == type)
            count[t->period * P + imin(R->X.pool->route[t->route].span, P - 1)]++;
    }
    /* the classes without voyage get a small objective too, so that the cut bears on them as much as possible */
    for (int i = 0; i < w->m; i++)
        weight[i] = count[w->cls[i]] + 1e-4;
    int status;
    GRB_CHECK(R->env, GRBsetdblattrarray(w->lp, GRB_DBL_ATTR_OBJ, 0, w->m, weight));
    GRB_CHECK(R->env, GRBoptimize(w->lp));
    GRB_CHECK(R->env, GRBgetintattr(w->lp, GRB_INT_ATTR_STATUS, &status));
    bool violated = false;
    if (status == GRB_OPTIMAL) {
        GRB_CHECK(R->env, GRBgetdblattrarray(w->lp, GRB_DBL_ATTR_X, 0, w->m, weight));
        /* exact validity whatever the solver's tolerances: no week may weigh more than 1 */
        double heaviest = 1.0, value = 0.0;
        for (int i = 0; i < w->m; i++)
            if (weight[i] < 1e-9)
                weight[i] = 0.0;
        for (int r = 0; r < w->rows; r++) {
            double sum = 0.0;
            for (int i = w->row_first[r]; i < w->row_first[r + 1]; i++)
                sum += weight[w->row_class[i]];
            heaviest = fmax(heaviest, sum);
        }
        for (int i = 0; i < w->m; i++) {
            weight[i] /= heaviest;
            value += weight[i] * count[w->cls[i]];
        }
        violated = value > plan->fleet[type] + 1e-6;
    }
    if (violated) {
        GrbRow *row = &R->row;
        row_clear(row);
        for (int i = 0; i < w->m; i++) {
            if (weight[i] == 0.0)
                continue;
            const int c = type * P * P + w->cls[i];
            for (int j = R->class_first[c]; j < R->class_first[c + 1]; j++)
                row_add(row, R->class_col[j], weight[i]);
        }
        row_add(row, R->y0 + type, -1.0);
        add_cut(R, cbdata, row, 0.0);
    }
    free(count);
    free(weight);
    return violated;
}

/* Cut (33): the voyages of type v marked as conflicting need @p needed vessels, and every one that
 * is missing saves at most one: y_v >= needed - (number of them missing). A voyage alone in its
 * class (day, span) is represented by c_vtl, so any voyage of the class stands for it; voyages
 * that share a class are represented by their own columns. (The paper uses c_vtl in both cases,
 * counted as many times as the class has voyages; a later plan with fewer voyages in the class
 * would then be cut off wrongly.) */
static void add_pattern_cut(Reduced *R, void *cbdata, const Plan *plan, int type, int needed)
{
    const int P = R->X.inst->num_periods;
    int *count = XCALLOC_ARRAY(int, P * P);
    for (int i = 0; i < plan->count; i++) {
        const Trip *t = &plan->trip[i];
        if (t->type == type && t->conflict)
            count[t->period * P + imin(R->X.pool->route[t->route].span, P - 1)]++;
    }
    GrbRow *row = &R->row;
    row_clear(row);
    int voyages = 0;
    for (int i = 0; i < plan->count; i++) {
        const Trip *t = &plan->trip[i];
        if (t->type != type || !t->conflict)
            continue;
        const int span = imin(R->X.pool->route[t->route].span, P - 1);
        row_add(row, count[t->period * P + span] == 1 ? c_var(R, type, t->period, span) : t->var, 1.0);
        voyages++;
    }
    row_add(row, R->y0 + type, -1.0);
    add_cut(R, cbdata, row, voyages - needed);
    free(count);
}

/* Adapters of the MIP models to ExactFallback. */
static bool model_relaxation(void *context, Plan *plan) { return model_berth_relaxation(context, plan); }

static int model_feasibility(void *context, Plan *plan, bool keep_assignment)
{
    return model_subproblem(context, plan, keep_assignment);
}

/* Writes the best plan to the solution file (through a temporary file, so that the file is always complete). */
static void save_best(Reduced *R)
{
    if (R->output == nullptr || R->best.count == 0)
        return;
    if (!sep_verify(&R->X, &R->best))
        die("internal error: an accepted plan fails the direct check of its schedule");
    char tmp[4096];
    if (snprintf(tmp, sizeof tmp, "%s.tmp", R->output) >= (int)sizeof tmp)
        die("solution file name too long: '%s'", R->output);
    FILE *f = fopen(tmp, "w");
    if (f == nullptr)
        die("cannot write '%s'", tmp);
    plan_print(&R->best, f);
    if (fclose(f) != 0 || rename(tmp, R->output) != 0)
        die("cannot write '%s'", R->output);
}

/* Prints the best plan's cost and the best bound so far, and flushes the output. */
static void print_progress(Reduced *R)
{
    const double t = wall_time() - R->wall0;
    if (R->best_obj < INFINITY)
        printf("progress: round %d, %.0f s, UB %.6f, LB %.6f, gap %.4f%%\n", R->round + 1, t, R->best_obj / COST_SCALE,
               R->lb / COST_SCALE, 100.0 * (R->best_obj - R->lb) / fmax(fabs(R->best_obj), 1e-9));
    else
        printf("progress: round %d, %.0f s, no plan yet, LB %.6f\n", R->round + 1, t, R->lb / COST_SCALE);
    fflush(stdout);
    R->next_progress = wall_time() + PROGRESS_SECONDS;
}

/* Solves the LP relaxation of the model, in place (no copy of the model): the variables are made continuous and
 * the constraints c_vtl = max(x) are removed, since a general constraint would make Gurobi solve a MIP. Before
 * the first cut the c_vtl appear in no other row and have no cost, so the LP value is that of the model. Presolve
 * is off, since on the largest models it alone can take hours; barrier without crossover solves these LPs (a
 * few hundred rows, millions of columns) quickly, and its value is lowered by ten times the barrier tolerance to
 * stay a valid bound. Returns the bound, or -INFINITY when the LP was not solved within @p limit seconds. The
 * model (types, constraints c_vtl = max(x)) and the parameters are restored. */
static double root_relaxation(Reduced *R, GRBmodel *m, double limit)
{
    const Instance *in = R->X.inst;
    const int P = in->num_periods, T = in->num_types;
    GRBenv *e = GRBgetenv(m);
    R->relaxing = true;
    int ngen;
    GRB_CHECK(R->env, GRBgetintattr(m, GRB_INT_ATTR_NUMGENCONSTRS, &ngen));
    if (ngen != T * P * P || R->num_cuts != 0)
        die("internal error: the root relaxation expects the model as built");
    int *gen = XMALLOC_ARRAY(int, ngen);
    for (int g = 0; g < ngen; g++)
        gen[g] = g;
    GRB_CHECK(R->env, GRBdelgenconstrs(m, ngen, gen));
    free(gen);
    char *vtype = XMALLOC_ARRAY(char, R->ncols), *cont = XMALLOC_ARRAY(char, R->ncols);
    GRB_CHECK(R->env, GRBgetcharattrarray(m, GRB_CHAR_ATTR_VTYPE, 0, R->ncols, vtype));
    memset(cont, GRB_CONTINUOUS, (size_t)R->ncols);
    GRB_CHECK(R->env, GRBsetcharattrarray(m, GRB_CHAR_ATTR_VTYPE, 0, R->ncols, cont));
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_PRESOLVE, 0));
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER));
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_CROSSOVER, 0));
    GRB_CHECK(R->env, GRBsetdblparam(e, GRB_DBL_PAR_TIMELIMIT, limit));
    GRB_CHECK(R->env, GRBoptimize(m));
    int status;
    double value = -INFINITY, tol;
    GRB_CHECK(R->env, GRBgetintattr(m, GRB_INT_ATTR_STATUS, &status));
    GRB_CHECK(R->env, GRBgetdblparam(e, GRB_DBL_PAR_BARCONVTOL, &tol));
    if (status == GRB_OPTIMAL) {
        GRB_CHECK(R->env, GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, &value));
        value -= 10 * tol * fmax(fabs(value), 1.0);
    }
    GRB_CHECK(R->env, GRBsetcharattrarray(m, GRB_CHAR_ATTR_VTYPE, 0, R->ncols, vtype));
    for (int k = 0; k < T; k++) /* as in build_model(): the class columns are its operands */
        for (int t = 0; t < P; t++)
            for (int l = 0; l < P; l++) {
                const int c = (k * P + t) * P + l;
                GRB_CHECK(R->env,
                          GRBaddgenconstrMax(m, nullptr, c_var(R, k, t, l), R->class_first[c + 1] - R->class_first[c],
                                             &R->class_col[R->class_first[c]], 0.0));
            }
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_PRESOLVE, -1));
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_METHOD, -1));
    GRB_CHECK(R->env, GRBsetintparam(e, GRB_INT_PAR_CROSSOVER, -1));
    GRB_CHECK(R->env, GRBupdatemodel(m));
    free(vtype);
    free(cont);
    R->relaxing = false;
    return value;
}

/* Separation procedure (Algorithm 1) for the integer solution in R->sol:
 * runs the checks of separation.h in order and adds the corresponding cut when one fails;
 * a solution passing all checks is kept in R->best if it is the best so far. */
static void check_incumbent(Reduced *R, void *cbdata)
{
    ExactContext *X = &R->X;
    const Instance *in = X->inst;
    Plan *plan = &R->plan;
    plan_clear(plan);
    for (int k = 0; k < in->num_types; k++)
        plan->fleet[k] = (int)lround(R->sol[R->y0 + k]);
    for (int j = R->nx - 1; j >= 0; --j)
        if (R->sol[j] > SELECTED)
            plan_add(plan, (Trip){.route = R->col[j].route,
                                  .period = R->col[j].period,
                                  .type = R->col[j].type,
                                  .vessel = -1,
                                  .berth = -1,
                                  .position = -1,
                                  .var = j});
    X->incumbents++;
    plan_big_m(plan);

    if (!sep_berth_relaxation(X, plan)) {
        add_nogood(R, cbdata, plan, true);
        X->berth_cuts++;
        return;
    }
    bool coloured = true;
    for (int k = in->num_types - 1; k >= 0; --k) {
        int needed;
        if (!sep_colouring(X, plan, k, &needed)) {
            if (add_week_cut(R, cbdata, plan, k)) {
                R->week_cuts++;
            } else {
                add_pattern_cut(R, cbdata, plan, k, needed);
                R->pattern_cuts++;
            }
            X->colouring_cuts++;
            coloured = false;
        }
    }
    if (!coloured)
        return;
    if (sep_subproblem(X, plan, true)) {
        /* the cost of the plan itself: the solver's value of a solution can be slightly off (its
         * variables are integral only within a tolerance) */
        const double cost = plan_cost(plan);
        if (cost < R->best_obj) {
            R->best_obj = cost;
            plan_copy(&R->best, plan);
            save_best(R);
        }
        return;
    }
    plan_copy(&R->reduced, plan);
    sep_reduce(X, &R->reduced);
    add_nogood(R, cbdata, &R->reduced, false);
    X->subproblem_cuts++;
}

/* Gurobi callback: calls check_incumbent() for every new integer solution (MIPSOL) and prints the progress
 * lines (MIP). Returns 0. */
static int callback(GRBmodel *model, void *cbdata, int where, void *usrdata)
{
    (void)model;
    Reduced *R = usrdata;
    if (R->relaxing)
        return 0;
    if (where == GRB_CB_MIP && R->progress && wall_time() >= R->next_progress) {
        double bound;
        if (GRBcbget(cbdata, where, GRB_CB_MIP_OBJBND, &bound) == 0)
            R->lb = fmax(R->lb, bound);
        print_progress(R);
        return 0;
    }
    if (where != GRB_CB_MIPSOL)
        return 0;
    GRB_CHECK(R->env, GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, R->sol));
    const double t0 = wall_time();
    check_incumbent(R, cbdata);
    R->X.check_seconds += wall_time() - t0;
    return 0;
}

/* ------------------------------------------------------------------------- */
/*  Driver                                                                   */
/* ------------------------------------------------------------------------- */

/* Adds the stored cuts to the model (when m is not null) and empties the store. */
static void flush_stored_cuts(Reduced *R, GRBmodel *m)
{
    for (int i = 0; i < R->num_cuts; i++) {
        StoredCut *c = &R->cuts[i];
        if (m)
            GRB_CHECK(R->env, GRBaddconstr(m, c->len, c->ind, c->val, GRB_LESS_EQUAL, c->rhs, nullptr));
        free(c->ind);
        free(c->val);
    }
    R->num_cuts = 0;
    if (m)
        GRB_CHECK(R->env, GRBupdatemodel(m));
}

int solve_reduced_model(const Instance *in, const ExactParams *par, bool quiet, const char *output)
{
    const double wall0 = wall_time(), cpu0 = process_cpu_time();
    if (!quiet)
        inst_print_summary(in);
    RoutePool *pool = route_pool_build(in, par->max_route_size);
    if (!quiet)
        route_pool_print_stats(pool);

    Reduced R = {.env = grb_open_env(par->solver_log),
                 .best_obj = INFINITY,
                 .output = output,
                 .progress = !quiet,
                 .wall0 = wall0,
                 .lb = -INFINITY};
    exact_context_init(&R.X, in, par, pool);
    R.week = XCALLOC_ARRAY(struct WeekModel, in->num_types);
    /* the models of the paper answer when a scheduling algorithm gives up */
    const int sub = par->subproblem_threads;
    ModelContext models = {.inst = in,
                           .pool = pool,
                           .env = R.env,
                           .threads = par->threads > 0 && (sub == 0 || sub > par->threads) ? par->threads : sub};
    rng_seed(&models.rng, 1, 0);
    R.X.fallback = (ExactFallback){model_relaxation, model_feasibility, &models};
    plan_init(&R.plan, in, pool);
    plan_init(&R.reduced, in, pool);
    plan_init(&R.best, in, pool);
    build_columns(&R);
    GRBmodel *m = build_model(&R);
    R.sol = XMALLOC_ARRAY(double, R.ncols);
    GRB_CHECK(R.env, GRBsetcallbackfunc(m, callback, &R));
    if (!quiet)
        printf("reduced model: %d voyage columns, %d columns in total%s\n", R.nx, R.ncols,
               installation_fixed_on_day_one(in) >= 0 ? " (one installation fixed on day 1 by symmetry)" : "");

    /* The bounds are kept across the rounds: the best accepted plan is R.best (the callback keeps it),
     * and since every cut is valid for the problem, the bound of every round is a valid lower bound.
     * (Adding the stored cuts after a round discards the solver's solution and bound.) */
    int status = GRB_LOADED;
    double gap = INFINITY, runtime = 0, lb = -INFINITY;
    GRB_CHECK(R.env, GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_MIPGAP, par->gap));
    /* a valid bound first, even if the run is stopped later (memory, job limit) */
    {
        const double t = wall_time();
        lb = R.lb = root_relaxation(&R, m, par->time_limit - (t - wall0));
        if (!quiet) {
            printf("root relaxation: LB %.6f, %.1f s\n", lb / COST_SCALE, wall_time() - t);
            fflush(stdout);
        }
    }
    for (int round = 0; !par->root_only; round++) {
        R.round = round;
        const int nlimits = (int)(sizeof ROUND_LIMITS / sizeof ROUND_LIMITS[0]);
        const double remaining = par->time_limit - (wall_time() - wall0);
        if (remaining <= 0)
            break;
        const double limit = fmin(ROUND_LIMITS[imin(round, nlimits - 1)], remaining);
        GRB_CHECK(R.env, GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, limit));
        GRB_CHECK(R.env, GRBoptimize(m));
        double t, bound;
        GRB_CHECK(R.env, GRBgetintattr(m, GRB_INT_ATTR_STATUS, &status));
        GRB_CHECK(R.env, GRBgetdblattr(m, GRB_DBL_ATTR_RUNTIME, &t));
        runtime += t;
        if (GRBgetdblattr(m, GRB_DBL_ATTR_OBJBOUND, &bound) == 0)
            lb = fmax(lb, bound);
        R.lb = fmax(R.lb, lb);
        if (R.best_obj < INFINITY)
            gap = (R.best_obj - lb) / fmax(fabs(R.best_obj), 1e-9);
        if (!quiet) {
            printf("round %d: status %d, %.1f s, gap %.4f%%, %d cuts to add\n", round + 1, status, t, 100 * gap,
                   R.num_cuts);
            fflush(stdout);
        }
        if (status == GRB_OPTIMAL || status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD || gap <= par->gap)
            break;
        if (status != GRB_TIME_LIMIT || par->time_limit - (wall_time() - wall0) < 1.0)
            break;
        /* restart with the cuts as ordinary constraints, starting from the incumbent */
        int sols = 0;
        GRB_CHECK(R.env, GRBgetintattr(m, GRB_INT_ATTR_SOLCOUNT, &sols));
        double *x = nullptr;
        if (sols > 0) {
            x = XMALLOC_ARRAY(double, R.ncols);
            GRB_CHECK(R.env, GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, R.ncols, x));
        }
        flush_stored_cuts(&R, m);
        if (x) {
            GRB_CHECK(R.env, GRBsetdblattrarray(m, GRB_DBL_ATTR_START, 0, R.ncols, x));
            free(x);
        }
    }
    if (par->root_only && !quiet)
        printf("stopped after the root relaxation (--root-only)\n");
    const double ub = R.best_obj;
    const bool solved = ub < INFINITY;
    const ExactContext *X = &R.X;
    /* optimal: the last round was solved, or the bound kept across the rounds is within the gap */
    const bool closed = status == GRB_OPTIMAL || (solved && gap <= par->gap);
    printf("branch-and-cut: %s, UB %.6f, LB %.6f, gap %.4f%%, solver %.1f s, wall %.1f s, cpu %.1f s\n",
           closed   ? "optimal"
           : solved ? "feasible"
                    : "no solution",
           ub / COST_SCALE, lb / COST_SCALE, solved ? 100.0 * (ub - lb) / fmax(fabs(ub), 1e-9) : INFINITY, runtime,
           wall_time() - wall0, process_cpu_time() - cpu0);
    printf("separation: %ld incumbents checked in %.1f s, cuts: %ld berth, %ld colouring (%ld weighted, %ld on "
           "patterns), %ld subproblem (%ld subproblems solved)\n",
           X->incumbents, X->check_seconds, X->berth_cuts, X->colouring_cuts, R.week_cuts, R.pattern_cuts,
           X->subproblem_cuts, X->subproblems);
    if (X->fallbacks > 0)
        printf("separation: %ld checks were answered by the MIP models (the scheduling algorithms reached their "
               "work limit)\n",
               X->fallbacks);
    if (X->unproven_subproblems > 0)
        printf("warning: %ld subproblems were left without answer and taken as infeasible; "
               "the result may not be optimal\n",
               X->unproven_subproblems);
    flush_stored_cuts(&R, nullptr);
    const bool found = R.best.count > 0;
    if (found && !sep_verify(&R.X, &R.best))
        die("internal error: the best plan fails the direct check of its schedule");
    if (found && !quiet) {
        printf("\nbest solution:\n");
        plan_print(&R.best, stdout);
    }
    if (found)
        save_best(&R);

    GRBfreemodel(m);
    for (int k = 0; k < in->num_types; k++) {
        if (R.week[k].lp)
            GRBfreemodel(R.week[k].lp);
        free(R.week[k].cls);
        free(R.week[k].row_first);
        free(R.week[k].row_class);
    }
    free(R.week);
    GRBfreeenv(R.env);
    exact_context_free(&R.X);
    row_free(&R.row);
    free(R.cuts);
    free(R.sol);
    free(R.col);
    free(R.z0);
    free(R.class_first);
    free(R.class_col);
    plan_free(&R.plan);
    plan_free(&R.reduced);
    plan_free(&R.best);
    route_pool_free(pool);
    return found ? EXIT_SUCCESS : EXIT_FAILURE;
}
