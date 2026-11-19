// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file complete_model.c
 * @brief Complete model (paper, Section 2.1) with lazy overlap cuts (see models.h).
 */
#include "exact/models.h"

#include "common/timer.h"
#include "common/util.h"
#include "exact/grb.h"
#include "exact/plan.h"
#include "exact/routes.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { K = MAX_DEPARTURES_PER_BERTH };

static const double SELECTED = 0.1;

/** @brief Column layout of the complete model and callback state. */
typedef struct Complete {
    const Instance *in;
    const ExactParams *par;
    const RoutePool *pool;
    GRBenv *env;
    long **routes; /* [type] voyages the type can perform */
    int *nroutes;  /* [type] */
    int *x0;       /* [type] first x column of the type */
    int twb0, tlb0, y0, *vessel0, z0, ncols;
    int *zi; /* [installation] first z column */
    double *sol;
    Plan plan, best;
    BerthSchedule bs;
    double best_obj;
    long incumbents, cuts;
    GrbRow row;
} Complete;

/* x(k, n, day i, berth j, position p, r-th voyage of type k) */
static int x_var(const Complete *C, int k, int n, int i, int j, int p, int r)
{
    const Instance *in = C->in;
    return C->x0[k] + ((((n * in->num_periods + i) * in->num_berths + j) * K + p) * C->nroutes[k] + r);
}
static int twb_var(const Complete *C, int i, int j) { return C->twb0 + i * C->in->num_berths + j; }
static int tlb_var(const Complete *C, int i, int j, int p) { return C->tlb0 + (i * C->in->num_berths + j) * K + p; }
/* Vessel block of type k: y(n) for every vessel, then (twv, tlv, t) for every (vessel, day). */
static int y_var(const Complete *C, int k, int n) { return C->vessel0[k] + n; }
/* Column of the waiting time t^WV of vessel @p n of type @p k departing on day @p i (t^LV and t follow it). */
static int twv_var(const Complete *C, int k, int n, int i)
{
    return C->vessel0[k] + C->in->type[k].max_vessels + (n * C->in->num_periods + i) * 3;
}
static int tlv_var(const Complete *C, int k, int n, int i) { return twv_var(C, k, n, i) + 1; }
static int t_var(const Complete *C, int k, int n, int i) { return twv_var(C, k, n, i) + 2; }

/* Appends @p count columns with the same objective coefficient, bounds and type to @p m. */
static void add_vars(GRBmodel *m, int count, double obj, double lb, double ub, char type)
{
    for (int i = 0; i < count; i++)
        GRB_CHECK(GRBgetenv(m), GRBaddvar(m, 0, nullptr, nullptr, obj, lb, ub, type, nullptr));
}

/* Builds the complete model (1)-(21) for the voyage pool of @p C and records the column layout in @p C.
 * Returns the model, ready to be optimised with the overlap-cut callback. */
static GRBmodel *build_model(Complete *C)
{
    const Instance *in = C->in;
    const RoutePool *pool = C->pool;
    const int P = in->num_periods, B = in->num_berths, T = in->num_types, N = in->num_installations;
    GRBmodel *m = grb_new_model(C->env, "complete", C->par->threads, C->par->solver_log);
    GRB_CHECK(C->env, GRBsetintparam(GRBgetenv(m), GRB_INT_PAR_LAZYCONSTRAINTS, 1));

    /* voyages per type */
    C->routes = XCALLOC_ARRAY(long *, T);
    C->nroutes = XCALLOC_ARRAY(int, T);
    C->x0 = XCALLOC_ARRAY(int, T);
    C->vessel0 = XCALLOC_ARRAY(int, T);
    int longest_load = 0;
    for (int k = 0; k < T; k++) {
        C->routes[k] = XMALLOC_ARRAY(long, imax((int)pool->count, 1));
        for (long r = 0; r < pool->count; r++)
            if (route_fits(&pool->route[r], k))
                C->routes[k][C->nroutes[k]++] = r;
    }
    for (long r = 0; r < pool->count; r++)
        longest_load = imax(longest_load, pool->route[r].load_time);
    /* big-M of the vessel waiting time for the second position: TOL - 1 plus the longest loading time */
    const double big_m = max_carry_days() + (double)longest_load / SECONDS_PER_DAY;

    /* columns */
    long next = 0;
    for (int k = 0; k < T; k++) {
        C->x0[k] = (int)next;
        const long n = (long)in->type[k].max_vessels * P * B * K * C->nroutes[k];
        if (next + n > INT32_MAX / 2)
            die("the complete model is too large for this instance (more than %d columns)", INT32_MAX / 2);
        for (int v = 0; v < in->type[k].max_vessels; v++)
            for (int i = 0; i < P; i++)
                for (int j = 0; j < B; j++)
                    for (int p = 0; p < K; p++)
                        for (int r = 0; r < C->nroutes[k]; r++) {
                            const double cost = (double)in->type[k].var_cost * pool->route[C->routes[k][r]].dist;
                            GRB_CHECK(C->env, GRBaddvar(m, 0, nullptr, nullptr, cost, 0.0, 1.0, GRB_BINARY, nullptr));
                        }
        next += n;
    }
    C->twb0 = (int)next;
    /* t^WB in [0, TOL - 1] (16); day P is the cyclic copy of day 0 */
    add_vars(m, (P + 1) * B, 0.0, 0.0, max_carry_days(), GRB_CONTINUOUS);
    next += (P + 1) * B;
    C->tlb0 = (int)next;
    add_vars(m, P * B * K, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS);
    next += P * B * K;
    for (int k = 0; k < T; k++) {
        C->vessel0[k] = (int)next;
        const int nv = in->type[k].max_vessels;
        add_vars(m, nv, in->type[k].fixed_cost, 0.0, 1.0, GRB_BINARY); /* y */
        for (int v = 0; v < nv; v++)
            for (int i = 0; i < P; i++) {
                add_vars(m, 2, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS); /* twv, tlv */
                add_vars(m, 1, 0.0, 0.0, P, GRB_INTEGER);               /* t */
            }
        next += nv + 3 * nv * P;
    }
    C->zi = XMALLOC_ARRAY(int, N);
    C->z0 = (int)next;
    for (int i = 0; i < N; i++) {
        C->zi[i] = (int)next;
        add_vars(m, in->inst[i].num_patterns, 0.0, 0.0, 1.0, GRB_BINARY);
        next += in->inst[i].num_patterns;
    }
    C->ncols = (int)next;
    GRB_CHECK(C->env, GRBupdatemodel(m));

    GrbRow row = {0};
    GrbRow *visit = XCALLOC_ARRAY(GrbRow, N * P);
    for (int k = 0; k < T; k++)
        for (int v = 0; v < in->type[k].max_vessels; v++)
            for (int i = 0; i < P; i++)
                for (int j = 0; j < B; j++)
                    for (int p = 0; p < K; p++)
                        for (int r = 0; r < C->nroutes[k]; r++) {
                            const long route = C->routes[k][r];
                            const int *s = route_stops(pool, route);
                            for (int a = 0; a < pool->route[route].size; a++)
                                row_add(&visit[s[a] * P + i], x_var(C, k, v, i, j, p, r), 1.0);
                        }
    for (int i = 0; i < N; i++) { /* (2) one pattern per installation */
        row_clear(&row);
        for (int p = 0; p < in->inst[i].num_patterns; p++)
            row_add(&row, C->zi[i] + p, 1.0);
        grb_add_row(m, &row, GRB_EQUAL, 1.0, nullptr);
    }
    for (int i = 0; i < N; i++) /* (3) visits follow the pattern */
        for (int t = 0; t < P; t++) {
            GrbRow *vr = &visit[i * P + t];
            for (int p = 0; p < in->inst[i].num_patterns; p++)
                if (inst_pattern_visits(&in->inst[i], P, p, t))
                    row_add(vr, C->zi[i] + p, -1.0);
            grb_add_row(m, vr, GRB_EQUAL, 0.0, nullptr);
            row_free(vr);
        }
    free(visit);

    for (int i = 0; i < P; i++) /* berth schedule (5)-(9) */
        for (int j = 0; j < B; j++) {
            for (int p = 0; p < K; p++) {
                row_clear(&row); /* (5) one voyage per departing position */
                for (int k = 0; k < T; k++)
                    for (int v = 0; v < in->type[k].max_vessels; v++)
                        for (int r = 0; r < C->nroutes[k]; r++)
                            row_add(&row, x_var(C, k, v, i, j, p, r), 1.0);
                grb_add_row(m, &row, GRB_LESS_EQUAL, 1.0, nullptr);
                row_clear(&row); /* (6) loading time of the position */
                row_add(&row, tlb_var(C, i, j, p), 1.0);
                for (int k = 0; k < T; k++)
                    for (int v = 0; v < in->type[k].max_vessels; v++)
                        for (int r = 0; r < C->nroutes[k]; r++)
                            row_add(&row, x_var(C, k, v, i, j, p, r),
                                    -(double)(float)((float)pool->route[C->routes[k][r]].load_time / SECONDS_PER_DAY));
                grb_add_row(m, &row, GRB_EQUAL, 0.0, nullptr);
            }
            row_clear(&row); /* (7) wait + loading <= TOL */
            row_add(&row, twb_var(C, i, j), 1.0);
            for (int p = 0; p < K; p++)
                row_add(&row, tlb_var(C, i, j, p), 1.0);
            grb_add_row(m, &row, GRB_LESS_EQUAL, tol_days(), nullptr);
            row_clear(&row); /* (8) carry-over to the next day */
            row_add(&row, twb_var(C, i + 1, j), 1.0);
            row_add(&row, twb_var(C, i, j), -1.0);
            for (int p = 0; p < K; p++)
                row_add(&row, tlb_var(C, i, j, p), -1.0);
            grb_add_row(m, &row, GRB_GREATER_EQUAL, -1.0, nullptr);
        }
    for (int j = 0; j < B; j++) { /* (9) cyclic horizon: day P is day 0 */
        row_clear(&row);
        row_add(&row, twb_var(C, P, j), 1.0);
        row_add(&row, twb_var(C, 0, j), -1.0);
        grb_add_row(m, &row, GRB_EQUAL, 0.0, nullptr);
    }

    for (int k = 0; k < T; k++)
        for (int v = 0; v < in->type[k].max_vessels; v++) {
            for (int i = 0; i < P; i++) {
                for (int j = 0; j < B; j++)
                    for (int p = 0; p < K; p++) { /* (10) waiting time of the vessel */
                        const double M = p == 0 ? max_carry_days() : big_m;
                        row_clear(&row);
                        row_add(&row, twv_var(C, k, v, i), 1.0);
                        row_add(&row, twb_var(C, i, j), -1.0);
                        for (int p2 = 0; p2 < p; p2++)
                            row_add(&row, tlb_var(C, i, j, p2), -1.0);
                        for (int r = 0; r < C->nroutes[k]; r++)
                            row_add(&row, x_var(C, k, v, i, j, p, r), -M);
                        grb_add_row(m, &row, GRB_GREATER_EQUAL, -M, nullptr);
                    }
                row_clear(&row); /* (11)-(12) loading and voyage times */
                row_add(&row, tlv_var(C, k, v, i), 1.0);
                for (int r = 0; r < C->nroutes[k]; r++) {
                    const double days = (double)(float)((float)pool->route[C->routes[k][r]].time / SECONDS_PER_DAY);
                    for (int j = 0; j < B; j++)
                        for (int p = 0; p < K; p++)
                            row_add(&row, x_var(C, k, v, i, j, p, r), -days);
                }
                grb_add_row(m, &row, GRB_EQUAL, 0.0, nullptr);
                row_clear(&row); /* (13) periods occupied by the vessel */
                row_add(&row, t_var(C, k, v, i), 1.0);
                row_add(&row, twv_var(C, k, v, i), -1.0);
                row_add(&row, tlv_var(C, k, v, i), -1.0);
                grb_add_row(m, &row, GRB_GREATER_EQUAL, 0.0, nullptr);
            }
            row_clear(&row); /* (14) within the horizon */
            for (int i = 0; i < P; i++)
                row_add(&row, t_var(C, k, v, i), 1.0);
            grb_add_row(m, &row, GRB_LESS_EQUAL, P, nullptr);
        }
    for (int i = 0; i < P; i++) /* (4) at most one active voyage per vessel and day */
        for (int k = 0; k < T; k++)
            for (int v = 0; v < in->type[k].max_vessels; v++) {
                row_clear(&row);
                for (int d = 0; d < P; d++)
                    for (int r = 0; r < C->nroutes[k]; r++)
                        if (route_active(in, &pool->route[C->routes[k][r]], d, i))
                            for (int j = 0; j < B; j++)
                                for (int p = 0; p < K; p++)
                                    row_add(&row, x_var(C, k, v, d, j, p, r), 1.0);
                row_add(&row, y_var(C, k, v), -1.0);
                grb_add_row(m, &row, GRB_LESS_EQUAL, 0.0, nullptr);
            }
    for (int k = 0; k < T; k++) /* vessels of a type are used in order (symmetry) */
        for (int v = in->type[k].max_vessels - 1; v > 0; --v) {
            row_clear(&row);
            row_add(&row, y_var(C, k, v), 1.0);
            row_add(&row, y_var(C, k, v - 1), -1.0);
            grb_add_row(m, &row, GRB_LESS_EQUAL, 0.0, nullptr);
        }
    row_free(&row);
    GRB_CHECK(C->env, GRBupdatemodel(m));
    return m;
}

/* Lazy cut (22) when some vessel cannot leave on time. */
static int callback(GRBmodel *model, void *cbdata, int where, void *usrdata)
{
    (void)model;
    Complete *C = usrdata;
    if (where != GRB_CB_MIPSOL)
        return 0;
    const Instance *in = C->in;
    const int P = in->num_periods, B = in->num_berths;
    GRB_CHECK(C->env, GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, C->sol));
    Plan *plan = &C->plan;
    plan_clear(plan);
    for (int k = 0; k < in->num_types; k++) {
        double used = 0;
        for (int v = 0; v < in->type[k].max_vessels; v++)
            used += C->sol[y_var(C, k, v)];
        plan->fleet[k] = (int)lround(used);
        for (int v = 0; v < in->type[k].max_vessels; v++)
            for (int i = 0; i < P; i++)
                for (int j = 0; j < B; j++)
                    for (int p = 0; p < K; p++)
                        for (int r = 0; r < C->nroutes[k]; r++) {
                            const int var = x_var(C, k, v, i, j, p, r);
                            if (C->sol[var] > SELECTED)
                                plan_add(plan, (Trip){.route = C->routes[k][r],
                                                      .period = (int8_t)i,
                                                      .type = (int8_t)k,
                                                      .vessel = (int8_t)v,
                                                      .berth = (int8_t)j,
                                                      .position = (int8_t)p,
                                                      .var = var});
                        }
    }
    C->incumbents++;
    if (plan_check_overlaps(plan, &C->bs)) {
        const double cost = plan_cost(plan); /* (the solver's value can be slightly off) */
        if (cost < C->best_obj) {
            C->best_obj = cost;
            plan_copy(&C->best, plan);
        }
        return 0;
    }
    row_clear(&C->row);
    for (int t = 0; t < plan->count; t++)
        if (plan->trip[t].conflict)
            row_add(&C->row, plan->trip[t].var, 1.0);
    GRB_CHECK(C->env, GRBcblazy(cbdata, C->row.len, C->row.ind, C->row.val, GRB_LESS_EQUAL, C->row.len - 1.0));
    C->cuts++;
    return 0;
}

int solve_complete_model(const Instance *in, const ExactParams *par, bool quiet, const char *output)
{
    const double wall0 = wall_time(), cpu0 = process_cpu_time();
    if (!quiet)
        inst_print_summary(in);
    RoutePool *pool = route_pool_build(in, par->max_route_size);
    if (!quiet)
        route_pool_print_stats(pool);
    Complete C = {.in = in, .par = par, .pool = pool, .env = grb_open_env(par->solver_log), .best_obj = INFINITY};
    plan_init(&C.plan, in, pool);
    plan_init(&C.best, in, pool);
    berth_schedule_init(&C.bs, in);
    GRBmodel *m = build_model(&C);
    C.sol = XMALLOC_ARRAY(double, C.ncols);
    if (!quiet)
        printf("complete model: %d columns\n", C.ncols);
    GRB_CHECK(C.env, GRBsetcallbackfunc(m, callback, &C));
    GRB_CHECK(C.env, GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, par->time_limit));
    GRB_CHECK(C.env, GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_MIPGAP, par->gap));
    GRB_CHECK(C.env, GRBoptimize(m));

    int status, sols;
    double ub = INFINITY, lb = -INFINITY, runtime;
    GRB_CHECK(C.env, GRBgetintattr(m, GRB_INT_ATTR_STATUS, &status));
    GRB_CHECK(C.env, GRBgetintattr(m, GRB_INT_ATTR_SOLCOUNT, &sols));
    GRB_CHECK(C.env, GRBgetdblattr(m, GRB_DBL_ATTR_RUNTIME, &runtime));
    if (sols > 0)
        ub = C.best_obj;
    if (GRBgetdblattr(m, GRB_DBL_ATTR_OBJBOUND, &lb) != 0)
        lb = -INFINITY;
    printf("complete model: %s, UB %.6f, LB %.6f, gap %.4f%%, solver %.1f s, wall %.1f s, cpu %.1f s\n",
           status == GRB_OPTIMAL ? "optimal"
           : sols > 0            ? "feasible"
                                 : "no solution",
           ub / COST_SCALE, lb / COST_SCALE, sols > 0 ? 100.0 * (ub - lb) / fmax(fabs(ub), 1e-9) : INFINITY, runtime,
           wall_time() - wall0, process_cpu_time() - cpu0);
    printf("overlap cuts: %ld (%ld incumbents checked)\n", C.cuts, C.incumbents);
    const bool found = C.best.count > 0;
    plan_compact_positions(&C.best);
    if (found && !quiet) {
        printf("\nbest solution:\n");
        plan_print(&C.best, stdout);
    }
    if (found && output) {
        FILE *f = fopen(output, "w");
        if (f == nullptr)
            die("cannot write '%s'", output);
        plan_print(&C.best, f);
        fclose(f);
    }
    GRBfreemodel(m);
    GRBfreeenv(C.env);
    for (int k = 0; k < in->num_types; k++)
        free(C.routes[k]);
    free(C.routes);
    free(C.nroutes);
    free(C.x0);
    free(C.vessel0);
    free(C.zi);
    free(C.sol);
    row_free(&C.row);
    plan_free(&C.plan);
    plan_free(&C.best);
    berth_schedule_free(&C.bs);
    route_pool_free(pool);
    return found ? EXIT_SUCCESS : EXIT_FAILURE;
}
