// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file separation.c
 * @brief Berth relaxation, fleet colouring, feasibility subproblem and Reduce
 *        (see separation.h).
 */
#include "exact/separation.h"

#include "common/util.h"

#include <stdlib.h>
#include <string.h>

void exact_context_init(ExactContext *X, const Instance *inst, const ExactParams *par, const RoutePool *pool)
{
    *X = (ExactContext){.inst = inst, .par = par, .pool = pool, .scheduler = scheduler_new()};
}

void exact_context_free(ExactContext *X)
{
    scheduler_free(X->scheduler);
    free(X->work);
    free(X->include);
    X->scheduler = nullptr;
    X->work = nullptr;
    X->include = nullptr;
    X->work_capacity = 0;
}

/** Reduce stops when the fallback had to answer this many of its subproblems (they can be slow). */
enum { REDUCE_MAX_FALLBACKS = 3 };

/* Scratch arrays of n entries each. */
enum {
    WORK_DAY,
    WORK_TYPE,
    WORK_LOAD,
    WORK_TIME,
    WORK_VESSEL,
    WORK_BERTH,
    WORK_POSITION,
    WORK_START,
    WORK_LENGTH,
    WORK_ARRAYS
};

static int *work(const ExactContext *X, int which) { return &X->work[which * X->work_capacity]; }

/* The plan as a scheduling problem (the arrays live in the context until the next call). */
static ScheduleProblem plan_problem(ExactContext *X, const Plan *plan)
{
    const Instance *in = X->inst;
    const int n = plan->count;
    if (n > X->work_capacity) {
        X->work_capacity = imax(2 * X->work_capacity, imax(n, 64));
        free(X->work);
        free(X->include);
        X->work = XMALLOC_ARRAY(int, WORK_ARRAYS * X->work_capacity);
        X->include = XMALLOC_ARRAY(bool, X->work_capacity);
    }
    int *day = work(X, WORK_DAY), *type = work(X, WORK_TYPE), *load = work(X, WORK_LOAD), *time = work(X, WORK_TIME);
    for (int r = 0; r < n; r++) {
        const Trip *t = &plan->trip[r];
        const Route *route = &plan->pool->route[t->route];
        day[r] = t->period;
        type[r] = t->type;
        load[r] = route->load_time;
        time[r] = route->time;
    }
    return (ScheduleProblem){.periods = in->num_periods,
                             .berths = in->num_berths,
                             .positions = MAX_DEPARTURES_PER_BERTH,
                             .day_length = SECONDS_PER_DAY,
                             .max_carry = BERTH_TOLERANCE - SECONDS_PER_DAY,
                             .num_types = in->num_types,
                             .fleet = plan->fleet,
                             .n = n,
                             .day = day,
                             .type = type,
                             .load = load,
                             .time = time};
}

/* ------------------------------------------------------------------------- */
/*  1. Berth relaxation                                                      */
/* ------------------------------------------------------------------------- */

bool sep_berth_relaxation(ExactContext *X, Plan *plan)
{
    const ScheduleProblem p = plan_problem(X, plan);
    const int n = plan->count;
    for (int r = 0; r < n; r++)
        plan->trip[r].conflict = false;
    if (schedule_berths(X->scheduler, &p, nullptr, nullptr, nullptr, nullptr))
        return true;
    if (scheduler_gave_up(X->scheduler)) {
        if (X->fallback.berth_relaxation == nullptr)
            return true; /* not settled here: the subproblem decides (and reports it) */
        X->fallbacks++;
        return X->fallback.berth_relaxation(X->fallback.context, plan);
    }
    /* The voyages of a few consecutive days cannot be scheduled on their own; among them, a minimal
     * infeasible subset: voyages are dropped one by one while the rest stays infeasible. */
    int first, days;
    scheduler_conflict_days(X->scheduler, &first, &days);
    bool *in = X->include;
    for (int r = 0; r < n; r++)
        in[r] = (p.day[r] - first + p.periods) % p.periods < days;
    for (int r = 0; r < n; r++) {
        if (!in[r])
            continue;
        in[r] = false;
        if (schedule_berths(X->scheduler, &p, in, nullptr, nullptr, nullptr) || scheduler_gave_up(X->scheduler))
            in[r] = true; /* needed (or not proven superfluous) */
    }
    for (int r = 0; r < n; r++)
        plan->trip[r].conflict = in[r];
    return false;
}

/* ------------------------------------------------------------------------- */
/*  2. Fleet colouring                                                       */
/* ------------------------------------------------------------------------- */

/* Can @p fleet vessels perform the voyages of @p type with include[r]? */
static bool fleet_suffices(ExactContext *X, const Plan *plan, int type, const bool *include, int fleet)
{
    int *start = work(X, WORK_START), *length = work(X, WORK_LENGTH), m = 0;
    for (int r = 0; r < plan->count; r++) {
        const Trip *t = &plan->trip[r];
        if (t->type != type || !include[r])
            continue;
        start[m] = t->period;
        length[m++] = plan->pool->route[t->route].span + 1;
    }
    return schedule_vessels(X->scheduler, X->inst->num_periods, m, start, length, fleet, nullptr);
}

bool sep_colouring(ExactContext *X, Plan *plan, int type, int *needed)
{
    const int n = plan->count, P = X->inst->num_periods;
    (void)plan_problem(X, plan); /* sizes the scratch arrays */
    bool *in = X->include;
    int count = 0;
    for (int r = 0; r < n; r++) {
        in[r] = plan->trip[r].type == type;
        count += in[r];
    }
    *needed = plan->fleet[type];
    if (count == 0 || fleet_suffices(X, plan, type, in, plan->fleet[type]))
        return true;
    for (int r = 0; r < n; r++)
        if (in[r])
            plan->trip[r].conflict = false;
    for (int r = 0; r < n; r++)
        if (in[r] && plan->pool->route[plan->trip[r].route].span + 1 > P) {
            /* a voyage longer than the week fits no vessel: it is the conflict on its own */
            plan->trip[r].conflict = true;
            *needed = plan->fleet[type] + 1;
            return false;
        }
    /* the smallest sufficient fleet (one vessel per voyage always works) */
    int fleet = plan->fleet[type] + 1;
    while (!fleet_suffices(X, plan, type, in, fleet))
        fleet++;
    /* a minimal set of voyages that needs it: drop voyages while one vessel fewer stays insufficient */
    for (int r = 0; r < n; r++) {
        if (!in[r])
            continue;
        in[r] = false;
        if (fleet_suffices(X, plan, type, in, fleet - 1))
            in[r] = true;
    }
    for (int r = 0; r < n; r++)
        if (plan->trip[r].type == type)
            plan->trip[r].conflict = in[r];
    *needed = fleet;
    return false;
}

/* ------------------------------------------------------------------------- */
/*  3. Feasibility subproblem                                                */
/* ------------------------------------------------------------------------- */

/* The subproblem: 1 feasible, 0 infeasible, -1 when neither the algorithm nor the fallback could tell. */
static int subproblem_answer(ExactContext *X, Plan *plan, bool keep_assignment)
{
    const Instance *in = X->inst;
    const int n = plan->count;
    X->subproblems++;
    if (n == 0)
        return 1;
    const ScheduleProblem p = plan_problem(X, plan);
    int *vessel = work(X, WORK_VESSEL), *berth = work(X, WORK_BERTH), *position = work(X, WORK_POSITION);
    const ScheduleStatus status = keep_assignment ? schedule_plan(X->scheduler, &p, vessel, berth, position)
                                                  : schedule_plan(X->scheduler, &p, nullptr, nullptr, nullptr);
    if (status == SCHEDULE_UNDECIDED) {
        if (X->fallback.subproblem == nullptr)
            return -1;
        X->fallbacks++;
        return X->fallback.subproblem(X->fallback.context, plan,
                                      keep_assignment); /* (it stores the assignment itself) */
    }
    if (status != SCHEDULE_FEASIBLE)
        return 0;
    if (keep_assignment) {
        for (int r = 0; r < n; r++) {
            plan->trip[r].vessel = (int8_t)vessel[r];
            plan->trip[r].berth = (int8_t)berth[r];
            plan->trip[r].position = (int8_t)position[r];
        }
        /* waiting times and slacks of the schedule, for the report */
        BerthSchedule bs;
        berth_schedule_init(&bs, in);
        plan_compact_positions(plan);
        for (int r = 0; r < n; r++) {
            const Trip *t = &plan->trip[r];
            *bs_load(&bs, t->period, t->berth, t->position) = trip_load_days(plan, t);
        }
        berth_schedule_waits(&bs);
        plan_fat_times(plan);
        plan_waits(plan, &bs);
        plan_lateness(plan);
        berth_schedule_free(&bs);
    }
    return 1;
}

bool sep_subproblem(ExactContext *X, Plan *plan, bool keep_assignment)
{
    const int answer = subproblem_answer(X, plan, keep_assignment);
    if (answer < 0)
        X->unproven_subproblems++;
    return answer > 0;
}

bool sep_verify(ExactContext *X, const Plan *plan)
{
    const ScheduleProblem p = plan_problem(X, plan);
    int *vessel = work(X, WORK_VESSEL), *berth = work(X, WORK_BERTH), *position = work(X, WORK_POSITION);
    for (int r = 0; r < plan->count; r++) {
        vessel[r] = plan->trip[r].vessel;
        berth[r] = plan->trip[r].berth;
        position[r] = plan->trip[r].position;
    }
    return schedule_check(&p, vessel, berth, position);
}

/* ------------------------------------------------------------------------- */
/*  Reduce (Algorithm 2)                                                     */
/* ------------------------------------------------------------------------- */

void sep_reduce(ExactContext *X, Plan *plan)
{
    const long fallbacks = X->fallbacks;
    for (int i = plan->count - 1; i >= 0 && X->fallbacks - fallbacks < REDUCE_MAX_FALLBACKS; --i) {
        /* move trip i to the end and drop it */
        const Trip t = plan->trip[i];
        plan->trip[i] = plan->trip[plan->count - 1];
        plan->trip[plan->count - 1] = t;
        plan->count--;
        if (subproblem_answer(X, plan, false) != 0)
            plan->count++; /* needed for the infeasibility (or not proven superfluous): put it back */
    }
}
