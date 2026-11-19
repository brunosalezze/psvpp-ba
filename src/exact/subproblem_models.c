// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file subproblem_models.c
 * @brief The MIP models of the berth relaxation and of the feasibility subproblem
 *        (see subproblem_models.h).
 */
#include "exact/subproblem_models.h"

#include "common/timer.h"
#include "common/util.h"

#include <stdlib.h>
#include <string.h>

enum { K = MAX_DEPARTURES_PER_BERTH };

/** Iteration limit of the AV/AB decomposition. */
enum { SUBPROBLEM_MAX_ITERATIONS = 10000 };
/** A cached conflict set is reused at most this many times. */
enum { TABU_MAX_REUSE = 50 };
/** Excess of the carried waiting time (days) below which a berth schedule counts as feasible. */
static const double BERTH_EXCESS_TOLERANCE = 0.0001;

/* ------------------------------------------------------------------------- */
/*  Berth assignment model (shared by the three berth checks)                */
/* ------------------------------------------------------------------------- */

typedef enum BerthVariant {
    BERTH_RELAXATION, /* minimise the largest excess wait; vessel timing ignored */
    BERTH_ASSIGNMENT, /* feasibility with the waiting time of every voyage bounded by its slack */
} BerthVariant;

typedef struct BerthModel {
    GRBmodel *m;
    int nvars;
    int tmw;  /* largest excess wait (relaxation only) */
    int *twb; /* [P*B] wait carried into (day, berth) */
    int *tlb; /* [P*B*K] loading time per (day, berth, position) */
    int *twv; /* [n] waiting time of each voyage (assignment only, -1 when excluded) */
    int *x;   /* [n*B*K] voyage -> (berth, position) (-1 when excluded) */
} BerthModel;

/* Releases a berth model and its index arrays. */
static void berth_model_free(BerthModel *bm)
{
    if (bm->m)
        GRBfreemodel(bm->m);
    free(bm->twb);
    free(bm->tlb);
    free(bm->twv);
    free(bm->x);
    *bm = (BerthModel){0};
}

/* Appends one column to the berth model and returns its index. */
static int add_var(BerthModel *bm, double obj, double lb, double ub, char type)
{
    GRB_CHECK(GRBgetenv(bm->m), GRBaddvar(bm->m, 0, nullptr, nullptr, obj, lb, ub, type, nullptr));
    return bm->nvars++;
}

/* Builds the berth model over the trips with include[r] (all when nullptr). */
static void berth_model_build(ModelContext *X, BerthModel *bm, const Plan *plan, const bool *include,
                              BerthVariant variant)
{
    const Instance *in = X->inst;
    const int P = in->num_periods, B = in->num_berths, n = plan->count;
    const bool assign = variant == BERTH_ASSIGNMENT;
    *bm = (BerthModel){.tmw = -1,
                       .twb = XMALLOC_ARRAY(int, P * B),
                       .tlb = XMALLOC_ARRAY(int, P * B * K),
                       .twv = XMALLOC_ARRAY(int, imax(n, 1)),
                       .x = XMALLOC_ARRAY(int, imax(n * B * K, 1))};
    bm->m = grb_new_model(X->env, assign ? "berth-assignment" : "berth-relaxation", X->threads, false);
    GRBenv *e = GRBgetenv(bm->m);
    if (assign) {
        GRB_CHECK(e, GRBsetintparam(e, GRB_INT_PAR_AGGFILL, 1000));
        GRB_CHECK(e, GRBsetintparam(e, GRB_INT_PAR_DEGENMOVES, 1));
        GRB_CHECK(e, GRBsetdblparam(e, GRB_DBL_PAR_HEURISTICS, 0));
        GRB_CHECK(e, GRBsetintparam(e, GRB_INT_PAR_STRONGCGCUTS, 2));
        GRB_CHECK(e, GRBsetintparam(e, GRB_INT_PAR_PRESOLVE, 2));
    }
#define IN(r) (include == nullptr || include[r])
    if (!assign)
        bm->tmw = add_var(bm, 1.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS);
    for (int i = 0; i < P; i++)
        for (int j = 0; j < B; j++) {
            bm->twb[i * B + j] = add_var(bm, 0.0, 0.0, assign ? max_carry_days() : GRB_INFINITY, GRB_CONTINUOUS);
            for (int k = 0; k < K; k++)
                bm->tlb[(i * B + j) * K + k] = add_var(bm, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS);
        }
    for (int r = 0; r < n; r++) {
        bm->twv[r] = -1;
        for (int jk = 0; jk < B * K; jk++)
            bm->x[r * B * K + jk] = -1;
        if (!IN(r))
            continue;
        if (assign)
            bm->twv[r] = add_var(bm, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS);
        for (int jk = 0; jk < B * K; jk++)
            bm->x[r * B * K + jk] = add_var(bm, 0.0, 0.0, 1.0, GRB_BINARY);
    }
    GRB_CHECK(e, GRBupdatemodel(bm->m));

    GrbRow row = {0};
    for (int i = 0; i < P; i++)
        for (int j = 0; j < B; j++) {
            if (!assign) { /* tmw >= twb - (TOL - 1) */
                row_clear(&row);
                row_add(&row, bm->tmw, 1.0);
                row_add(&row, bm->twb[i * B + j], -1.0);
                grb_add_row(bm->m, &row, GRB_GREATER_EQUAL, -max_carry_days(), nullptr);
            }
            for (int k = 0; k < K; k++) {
                row_clear(&row); /* at most one voyage per departing position */
                for (int r = n - 1; r >= 0; --r)
                    if (IN(r) && plan->trip[r].period == i)
                        row_add(&row, bm->x[(r * B + j) * K + k], 1.0);
                grb_add_row(bm->m, &row, GRB_LESS_EQUAL, 1.0, nullptr);
                row_clear(&row); /* loading time of the position */
                row_add(&row, bm->tlb[(i * B + j) * K + k], 1.0);
                for (int r = n - 1; r >= 0; --r)
                    if (IN(r) && plan->trip[r].period == i)
                        row_add(&row, bm->x[(r * B + j) * K + k], -trip_load_days(plan, &plan->trip[r]));
                grb_add_row(bm->m, &row, GRB_EQUAL, 0.0, nullptr);
            }
            row_clear(&row); /* carry-over: twb(t+1) >= twb(t) + sum_k tlb(t,k) - 1 */
            row_add(&row, bm->twb[((i + 1) % P) * B + j], 1.0);
            row_add(&row, bm->twb[i * B + j], -1.0);
            for (int k = 0; k < K; k++)
                row_add(&row, bm->tlb[(i * B + j) * K + k], -1.0);
            grb_add_row(bm->m, &row, GRB_GREATER_EQUAL, -1.0, nullptr);
        }
    for (int r = n - 1; r >= 0; --r) {
        if (!IN(r))
            continue;
        row_clear(&row); /* every voyage gets one berth and position */
        for (int jk = 0; jk < B * K; jk++)
            row_add(&row, bm->x[r * B * K + jk], 1.0);
        grb_add_row(bm->m, &row, GRB_EQUAL, 1.0, nullptr);
        if (!assign)
            continue;
        const int p = plan->trip[r].period;
        for (int j = 0; j < B; j++)
            for (int k = 0; k < K; k++) {
                /* twv >= twb + loads ahead in the queue, when the voyage takes (j, k) */
                const double big_m = k == 0 ? max_carry_days() : plan->trip[r].big_m;
                row_clear(&row);
                row_add(&row, bm->twv[r], 1.0);
                row_add(&row, bm->twb[p * B + j], -1.0);
                for (int k2 = 0; k2 < k; k2++)
                    row_add(&row, bm->tlb[(p * B + j) * K + k2], -1.0);
                row_add(&row, bm->x[(r * B + j) * K + k], -big_m);
                grb_add_row(bm->m, &row, GRB_GREATER_EQUAL, -big_m, nullptr);
            }
    }
    if (assign) { /* berths are identical: the first voyage of the earliest busy day uses berth 1 */
        int first = -1;
        for (int i = 0; i < P && first < 0; i++)
            for (int r = n - 1; r >= 0; --r)
                if (IN(r) && plan->trip[r].period == i) {
                    first = r;
                    break;
                }
        if (first >= 0) {
            row_clear(&row);
            for (int k = 0; k < K; k++)
                row_add(&row, bm->x[(first * B + 0) * K + k], 1.0);
            grb_add_row(bm->m, &row, GRB_EQUAL, 1.0, nullptr);
        }
    }
#undef IN
    row_free(&row);
    GRB_CHECK(e, GRBupdatemodel(bm->m));
}

/* Bounds the waiting time of every included voyage by its slack. */
static void berth_model_add_slacks(BerthModel *bm, const Plan *plan)
{
    for (int r = 0; r < plan->count; r++)
        if (bm->twv[r] >= 0) {
            int ind = bm->twv[r];
            double one = 1.0;
            GRB_CHECK(GRBgetenv(bm->m),
                      GRBaddconstr(bm->m, 1, &ind, &one, GRB_LESS_EQUAL, plan->trip[r].fat_time, nullptr));
        }
}

/* Optimises @p m and returns Gurobi's optimisation status (GRB_OPTIMAL, GRB_INFEASIBLE, ...). */
static int optimize_status(GRBmodel *m)
{
    int status;
    GRB_CHECK(GRBgetenv(m), GRBoptimize(m));
    GRB_CHECK(GRBgetenv(m), GRBgetintattr(m, GRB_INT_ATTR_STATUS, &status));
    return status;
}

/* ------------------------------------------------------------------------- */
/*  Berth relaxation                                                         */
/* ------------------------------------------------------------------------- */

/* Smallest achievable excess of the carried waiting time over TOL - 1 for the trips with
 * include[r] (all when nullptr), or a negative value when the model cannot be solved. When
 * @p sol is not null it receives the solution (bm->nvars values). */
static double berth_relaxation_excess(ModelContext *X, BerthModel *bm, const Plan *plan, const bool *include,
                                      double **sol)
{
    berth_model_build(X, bm, plan, include, BERTH_RELAXATION);
    if (optimize_status(bm->m) != GRB_OPTIMAL)
        return -1.0;
    double *x = XMALLOC_ARRAY(double, bm->nvars);
    GRB_CHECK(X->env, GRBgetdblattrarray(bm->m, GRB_DBL_ATTR_X, 0, bm->nvars, x));
    const double excess = x[bm->tmw];
    if (sol)
        *sol = x;
    else
        free(x);
    return excess;
}

bool model_berth_relaxation(ModelContext *X, Plan *plan)
{
    const Instance *in = X->inst;
    const int P = in->num_periods, B = in->num_berths, n = plan->count;
    BerthModel bm;
    double *sol = nullptr;
    const double excess = berth_relaxation_excess(X, &bm, plan, nullptr, &sol);
    for (int r = 0; r < n; r++)
        plan->trip[r].conflict = excess < 0; /* unsolved model: blame every voyage */
    if (excess >= 0 && excess <= BERTH_EXCESS_TOLERANCE) {
        free(sol);
        berth_model_free(&bm);
        return true;
    }
    if (excess > 0) {
        /* most occupied (day, berth) */
        double best = 0;
        int day = 0;
        for (int i = P - 1; i >= 0; --i)
            for (int j = B - 1; j >= 0; --j) {
                double use = sol[bm.twb[i * B + j]];
                for (int k = K - 1; k >= 0; --k)
                    use += sol[bm.tlb[(i * B + j) * K + k]];
                if (use > best) {
                    best = use;
                    day = i;
                }
            }
        /* that day and the preceding days that carried waiting time into it */
        for (int c = 0; c < P; c++, day = day == 0 ? P - 1 : day - 1) {
            for (int r = n - 1; r >= 0; --r)
                if (plan->trip[r].period == day)
                    plan->trip[r].conflict = true;
            bool carried = false;
            for (int j = B - 1; j >= 0; --j)
                carried |= sol[bm.twb[day * B + j]] > 0.001;
            if (!carried)
                break;
        }
    }
    free(sol);
    berth_model_free(&bm);

    /* The voyages of those days are only a guess of the cause: the solver may have traded waiting
     * time between days. The cut built from them is valid only if they cannot be scheduled on
     * their own, so check it, and blame every voyage of the plan otherwise. */
    bool *subset = XMALLOC_ARRAY(bool, imax(n, 1));
    bool all = true;
    for (int r = 0; r < n; r++)
        all &= subset[r] = plan->trip[r].conflict;
    if (!all) {
        const double alone = berth_relaxation_excess(X, &bm, plan, subset, nullptr);
        berth_model_free(&bm);
        if (alone >= 0 && alone <= BERTH_EXCESS_TOLERANCE)
            for (int r = 0; r < n; r++)
                plan->trip[r].conflict = true;
    }
    free(subset);
    return false;
}

/* ------------------------------------------------------------------------- */
/*  Vessel assignment model                                                  */
/* ------------------------------------------------------------------------- */

/* Builds the vessel assignment model AV over the trips of the given type
 * (all types when type < 0). x[r * max_fleet + n] is the column of trip r on vessel n. */
static GRBmodel *vessel_model_build(ModelContext *X, const Plan *plan, int type, int *x, int max_fleet, int threads)
{
    const Instance *in = X->inst;
    const int n = plan->count;
    GRBmodel *m = grb_new_model(X->env, "vessel-assignment", threads, false);
    int nv = 0;
    for (int r = n - 1; r >= 0; --r) {
        const Trip *t = &plan->trip[r];
        for (int v = 0; v < max_fleet; v++)
            x[r * max_fleet + v] = -1;
        if (type >= 0 && t->type != type)
            continue;
        for (int v = plan->fleet[t->type] - 1; v >= 0; --v) {
            GRB_CHECK(X->env, GRBaddvar(m, 0, nullptr, nullptr, 0.0, 0.0, 1.0, GRB_BINARY, nullptr));
            x[r * max_fleet + v] = nv++;
        }
    }
    GRB_CHECK(X->env, GRBupdatemodel(m));
    GrbRow row = {0};
    for (int r = n - 1; r >= 0; --r) { /* every voyage gets one vessel */
        const Trip *t = &plan->trip[r];
        if (type >= 0 && t->type != type)
            continue;
        row_clear(&row);
        for (int v = plan->fleet[t->type] - 1; v >= 0; --v)
            row_add(&row, x[r * max_fleet + v], 1.0);
        grb_add_row(m, &row, GRB_EQUAL, 1.0, nullptr);
    }
    for (int i = in->num_periods - 1; i >= 0; --i) /* a vessel performs one voyage at a time */
        for (int k = in->num_types - 1; k >= 0; --k) {
            if (type >= 0 && k != type)
                continue;
            for (int v = plan->fleet[k] - 1; v >= 0; --v) {
                row_clear(&row);
                for (int r = n - 1; r >= 0; --r) {
                    const Trip *t = &plan->trip[r];
                    if (t->type == k && route_active(in, &X->pool->route[t->route], t->period, i))
                        row_add(&row, x[r * max_fleet + v], 1.0);
                }
                if (row.len > 0)
                    grb_add_row(m, &row, GRB_LESS_EQUAL, 1.0, nullptr);
            }
        }
    row_free(&row);
    return m;
}

/* Largest fleet of any vessel type in the plan (at least 1): the stride of the vessel-assignment columns. */
static int max_fleet_of(const Plan *plan)
{
    int f = 1;
    for (int k = 0; k < plan->inst->num_types; k++)
        f = imax(f, plan->fleet[k]);
    return f;
}

/* ------------------------------------------------------------------------- */
/*  Feasibility subproblem                                                   */
/* ------------------------------------------------------------------------- */

/* A cached infeasible configuration: the conflicting trips and their slacks. */
typedef struct TabuEntry {
    float *fat_time;
    bool *conflict;
    int reuse;
} TabuEntry;

typedef struct TabuList {
    TabuEntry *entry;
    int count, capacity, n;
} TabuList;

/* Releases the entries of a tabu list. */
static void tabu_free(TabuList *t)
{
    for (int i = 0; i < t->count; i++) {
        free(t->entry[i].fat_time);
        free(t->entry[i].conflict);
    }
    free(t->entry);
}

/* Stores the current slacks and conflict flags of the plan's trips as a new tabu entry. */
static void tabu_add(TabuList *t, const Plan *plan)
{
    if (t->count == t->capacity) {
        t->capacity = t->capacity ? 2 * t->capacity : 16;
        t->entry = xrealloc(t->entry, sizeof(TabuEntry) * (size_t)t->capacity);
    }
    TabuEntry *e = &t->entry[t->count++];
    e->fat_time = XMALLOC_ARRAY(float, t->n);
    e->conflict = XMALLOC_ARRAY(bool, t->n);
    e->reuse = 0;
    for (int r = 0; r < t->n; r++) {
        e->fat_time[r] = plan->trip[r].fat_time;
        e->conflict[r] = plan->trip[r].conflict;
    }
}

/* A cached conflict set applies when no conflicting trip has more slack now
 * than when the set was found infeasible. Sets Trip::conflict on a hit. */
static bool tabu_lookup(TabuList *t, Plan *plan)
{
    for (int i = t->count - 1; i >= 0; --i) {
        TabuEntry *e = &t->entry[i];
        if (e->reuse >= TABU_MAX_REUSE)
            continue;
        bool match = true;
        for (int r = t->n - 1; r >= 0 && match; --r)
            if (e->conflict[r] && e->fat_time[r] < plan->trip[r].fat_time - 0.001f)
                match = false;
        if (match) {
            e->reuse++;
            for (int r = t->n - 1; r >= 0; --r)
                plan->trip[r].conflict = e->conflict[r];
            return true;
        }
    }
    return false;
}

/* Berth assignment over the conflicting trips only, with their slacks. */
static bool conflict_set_feasible(ModelContext *X, const Plan *plan, const bool *subset)
{
    BerthModel bm;
    berth_model_build(X, &bm, plan, subset, BERTH_ASSIGNMENT);
    berth_model_add_slacks(&bm, plan);
    const int status = optimize_status(bm.m);
    berth_model_free(&bm);
    return status == GRB_OPTIMAL;
}

/* Shrinks the berth-infeasible set of trips (with the current vessel
 * assignment) to a minimal infeasible subset, examined in random order. */
static void minimal_berth_conflict(ModelContext *X, Plan *plan)
{
    const int n = plan->count;
    bool *subset = XMALLOC_ARRAY(bool, n);
    int *order = XMALLOC_ARRAY(int, n);
    for (int r = 0; r < n; r++) {
        subset[r] = true;
        order[r] = r;
    }
    rng_shuffle_int(&X->rng, order, n);
    for (int i = n - 1; i >= 0; --i) {
        const int r = order[i];
        subset[r] = false;
        if (conflict_set_feasible(X, plan, subset))
            subset[r] = true;
    }
    for (int r = 0; r < n; r++)
        plan->trip[r].conflict = subset[r];
    free(subset);
    free(order);
}

/* The trip whose departure ends the time window of trip @p r: the next departure (cyclically) of
 * the vessel assigned to @p r. Returns -1 when the vessel performs no other trip. */
static int next_trip_of_vessel(const Plan *plan, int r)
{
    const Trip *t = &plan->trip[r];
    const int P = plan->inst->num_periods;
    int next = -1, best_gap = P + 1;
    for (int q = 0; q < plan->count; q++) {
        const Trip *u = &plan->trip[q];
        if (q == r || u->type != t->type || u->vessel != t->vessel)
            continue;
        const int gap = (u->period - t->period + P) % P;
        if (gap > 0 && gap < best_gap) {
            best_gap = gap;
            next = q;
        }
    }
    return next;
}

int model_subproblem(ModelContext *X, Plan *plan, bool keep_assignment)
{
    const Instance *in = X->inst;
    const int P = in->num_periods, B = in->num_berths, T = in->num_types, n = plan->count;
    if (n == 0)
        return 1;
    plan_big_m(plan);
    bool unproven = false;

    /* Trips active on each day per type; the busiest day fixes the vessel
     * numbering (symmetry breaking). */
    int *active = XCALLOC_ARRAY(int, T * P);
    int *busy_trips = XMALLOC_ARRAY(int, T * P * n);
    for (int r = n - 1; r >= 0; --r) {
        const Trip *t = &plan->trip[r];
        for (int d = 0; d < P; d++)
            if (route_active(in, &X->pool->route[t->route], t->period, d))
                busy_trips[(t->type * P + d) * n + active[t->type * P + d]++] = r;
    }

    const int mf = max_fleet_of(plan);
    int *xv = XMALLOC_ARRAY(int, n * mf);
    GRBmodel *av = vessel_model_build(X, plan, -1, xv, mf, X->threads);
    bool possible = true;
    GrbRow row = {0};
    for (int k = T - 1; k >= 0; --k) {
        int busiest = 0;
        for (int d = P - 1; d > 0; --d)
            if (active[k * P + d] > active[k * P + busiest])
                busiest = d;
        if (active[k * P + busiest] > plan->fleet[k]) {
            possible = false; /* more simultaneous voyages than vessels */
            break;
        }
        for (int a = active[k * P + busiest] - 1; a >= 0; --a) {
            const int r = busy_trips[(k * P + busiest) * n + a];
            row_clear(&row);
            row_add(&row, xv[r * mf + a], 1.0);
            grb_add_row(av, &row, GRB_EQUAL, 1.0, nullptr);
        }
    }
    free(active);
    free(busy_trips);

    BerthModel ab = {0};
    TabuList tabu = {.n = n};
    double *sol = XMALLOC_ARRAY(double, imax(n * mf, 1));
    bool *in_cut = XMALLOC_ARRAY(bool, n);
    bool feasible = false;
    for (int it = 0; possible && it < SUBPROBLEM_MAX_ITERATIONS; it++) {
        if (optimize_status(av) != GRB_OPTIMAL)
            break; /* no vessel assignment left */
        int nav;
        GRB_CHECK(X->env, GRBgetintattr(av, GRB_INT_ATTR_NUMVARS, &nav));
        GRB_CHECK(X->env, GRBgetdblattrarray(av, GRB_DBL_ATTR_X, 0, nav, sol));
        for (int r = 0; r < n; r++)
            for (int v = 0; v < plan->fleet[plan->trip[r].type]; v++)
                if (sol[xv[r * mf + v]] > 0.5)
                    plan->trip[r].vessel = (int8_t)v;
        plan_fat_times(plan);

        if (!tabu_lookup(&tabu, plan)) {
            berth_model_build(X, &ab, plan, nullptr, BERTH_ASSIGNMENT);
            berth_model_add_slacks(&ab, plan);
            if (optimize_status(ab.m) == GRB_OPTIMAL) {
                feasible = true;
                break;
            }
            berth_model_free(&ab);
            minimal_berth_conflict(X, plan);
            tabu_add(&tabu, plan);
        }
        /* Forbid this vessel assignment of the conflicting trips. Their berth conflict was
         * proved with the slack that each of them has now, which lasts until the next departure
         * of its vessel; that next voyage must therefore stay on the vessel for the cut to hold
         * (with it elsewhere, the conflicting trip could have more slack and fit). */
        for (int r = 0; r < n; r++)
            in_cut[r] = plan->trip[r].conflict;
        for (int r = 0; r < n; r++) {
            const int next = plan->trip[r].conflict ? next_trip_of_vessel(plan, r) : -1;
            if (next >= 0)
                in_cut[next] = true;
        }
        row_clear(&row);
        for (int r = n - 1; r >= 0; --r)
            if (in_cut[r])
                row_add(&row, xv[r * mf + plan->trip[r].vessel], 1.0);
        grb_add_row(av, &row, GRB_LESS_EQUAL, row.len - 1.0, nullptr);
        if (it + 1 == SUBPROBLEM_MAX_ITERATIONS)
            unproven = true;
    }

    if (feasible && keep_assignment) {
        double *x = XMALLOC_ARRAY(double, ab.nvars);
        GRB_CHECK(X->env, GRBgetdblattrarray(ab.m, GRB_DBL_ATTR_X, 0, ab.nvars, x));
        BerthSchedule bs;
        berth_schedule_init(&bs, in);
        for (int r = 0; r < n; r++)
            for (int jk = 0; jk < B * K; jk++)
                if (x[ab.x[r * B * K + jk]] > 0.5) {
                    plan->trip[r].berth = (int8_t)(jk / K);
                    plan->trip[r].position = (int8_t)(jk % K);
                }
        plan_compact_positions(plan);
        for (int r = 0; r < n; r++) {
            const Trip *t = &plan->trip[r];
            *bs_load(&bs, t->period, t->berth, t->position) = trip_load_days(plan, t);
        }
        berth_schedule_waits(&bs);
        plan_waits(plan, &bs);
        plan_lateness(plan);
        berth_schedule_free(&bs);
        free(x);
    }
    berth_model_free(&ab);
    tabu_free(&tabu);
    row_free(&row);
    free(sol);
    free(in_cut);
    free(xv);
    GRBfreemodel(av);
    return feasible ? 1 : unproven ? -1 : 0;
}
