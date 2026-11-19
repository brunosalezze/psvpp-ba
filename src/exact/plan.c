// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file plan.c
 * @brief Plans of trips and schedule computations (see plan.h).
 */
#include "exact/plan.h"

#include "common/util.h"

#include <stdlib.h>
#include <string.h>

/** Lateness below this value (days) is rounding noise. */
static const float LATE_TOLERANCE = 0.001f;

void plan_init(Plan *plan, const Instance *inst, const RoutePool *pool)
{
    *plan = (Plan){.inst = inst, .pool = pool, .fleet = XCALLOC_ARRAY(int, inst->num_types)};
}

void plan_free(Plan *plan)
{
    free(plan->trip);
    free(plan->fleet);
    *plan = (Plan){0};
}

void plan_clear(Plan *plan)
{
    plan->count = 0;
    memset(plan->fleet, 0, sizeof(int) * (size_t)plan->inst->num_types);
}

void plan_add(Plan *plan, Trip trip)
{
    if (plan->count == plan->capacity) {
        plan->capacity = plan->capacity ? 2 * plan->capacity : 32;
        plan->trip = xrealloc(plan->trip, sizeof(Trip) * (size_t)plan->capacity);
    }
    plan->trip[plan->count++] = trip;
}

void plan_copy(Plan *dst, const Plan *src)
{
    dst->count = 0;
    for (int i = 0; i < src->count; i++)
        plan_add(dst, src->trip[i]);
    memcpy(dst->fleet, src->fleet, sizeof(int) * (size_t)src->inst->num_types);
}

double plan_cost(const Plan *plan)
{
    const Instance *in = plan->inst;
    double cost = 0;
    for (int k = 0; k < in->num_types; k++)
        cost += (double)plan->fleet[k] * in->type[k].fixed_cost;
    for (int i = 0; i < plan->count; i++) {
        const Trip *t = &plan->trip[i];
        cost += (double)in->type[t->type].var_cost * plan->pool->route[t->route].dist;
    }
    return cost;
}

float trip_load_days(const Plan *plan, const Trip *t)
{
    return (float)plan->pool->route[t->route].load_time / SECONDS_PER_DAY;
}

void berth_schedule_init(BerthSchedule *bs, const Instance *inst)
{
    const int n = inst->num_periods * inst->num_berths;
    *bs = (BerthSchedule){
        .inst = inst, .wait = XCALLOC_ARRAY(float, n), .load = XCALLOC_ARRAY(float, n * MAX_DEPARTURES_PER_BERTH)};
}

void berth_schedule_free(BerthSchedule *bs)
{
    free(bs->wait);
    free(bs->load);
    *bs = (BerthSchedule){0};
}

/* Occupation of a berth on a day: carried wait plus the loading times. */
static float occupation(const BerthSchedule *bs, int p, int b)
{
    float sum = *bs_wait(bs, p, b);
    for (int k = MAX_DEPARTURES_PER_BERTH - 1; k >= 0; --k)
        sum += *bs_load(bs, p, b, k);
    return sum;
}

void berth_schedule_waits(BerthSchedule *bs)
{
    const Instance *in = bs->inst;
    const int P = in->num_periods;
    for (int p = 0; p < P; p++)
        for (int b = 0; b < in->num_berths; b++)
            *bs_wait(bs, p, b) = 0;
    /* Three cyclic sweeps propagate any carry-over around the week. */
    for (int sweep = 0; sweep < 3; sweep++)
        for (int b = in->num_berths - 1; b >= 0; --b) {
            for (int p = 0; p < P - 1; p++) {
                const float excess = occupation(bs, p, b) - 1;
                *bs_wait(bs, p + 1, b) = excess > 0 ? excess : 0;
            }
            const float excess = occupation(bs, P - 1, b) - 1;
            *bs_wait(bs, 0, b) = excess > 0 ? excess : 0;
        }
}

void plan_fat_times(Plan *plan)
{
    const Instance *in = plan->inst;
    const int P = in->num_periods;
    /* busy[type][vessel][day] = the vessel departs on that day */
    int max_vessels = 1;
    for (int k = 0; k < in->num_types; k++)
        max_vessels = imax(max_vessels, plan->fleet[k]);
    for (int i = 0; i < plan->count; i++)
        max_vessels = imax(max_vessels, plan->trip[i].vessel + 1);
    bool *busy = XCALLOC_ARRAY(bool, in->num_types * max_vessels * P);
#define BUSY(t, v, p) busy[((t) * max_vessels + (v)) * P + (p)]
    for (int i = 0; i < plan->count; i++)
        BUSY(plan->trip[i].type, plan->trip[i].vessel, plan->trip[i].period) = true;
    for (int i = 0; i < plan->count; i++) {
        Trip *t = &plan->trip[i];
        int gap = 1;
        for (int d = t->period + 1; gap <= P; d++, gap++) {
            if (d == P)
                d = 0;
            if (BUSY(t->type, t->vessel, d))
                break;
        }
        t->fat_time = (float)(gap * SECONDS_PER_DAY - plan->pool->route[t->route].time) / SECONDS_PER_DAY;
    }
#undef BUSY
    free(busy);
}

void plan_waits(Plan *plan, const BerthSchedule *bs)
{
    for (int i = 0; i < plan->count; i++) {
        Trip *t = &plan->trip[i];
        float w = *bs_wait(bs, t->period, t->berth);
        for (int k = t->position - 1; k >= 0; --k)
            w += *bs_load(bs, t->period, t->berth, k);
        t->wait = w;
    }
}

void plan_lateness(Plan *plan)
{
    for (int i = 0; i < plan->count; i++) {
        Trip *t = &plan->trip[i];
        t->lateness = t->wait > t->fat_time ? t->wait - t->fat_time : 0;
    }
}

void plan_big_m(Plan *plan)
{
    for (int i = 0; i < plan->count; i++) {
        int longest = 0; /* its own loading time counts: see plan.h */
        for (int j = 0; j < plan->count; j++)
            if (plan->trip[j].period == plan->trip[i].period)
                longest = imax(longest, plan->pool->route[plan->trip[j].route].load_time);
        plan->trip[i].big_m = (float)max_carry_days() + (float)longest / SECONDS_PER_DAY;
    }
}

bool plan_check_overlaps(Plan *plan, BerthSchedule *bs)
{
    const Instance *in = plan->inst;
    const int P = in->num_periods;
    for (int p = 0; p < P; p++)
        for (int b = 0; b < in->num_berths; b++)
            for (int k = 0; k < MAX_DEPARTURES_PER_BERTH; k++)
                *bs_load(bs, p, b, k) = 0;
    for (int i = plan->count - 1; i >= 0; --i) {
        Trip *t = &plan->trip[i];
        t->conflict = false;
        *bs_load(bs, t->period, t->berth, t->position) = trip_load_days(plan, t);
    }
    plan_fat_times(plan);
    berth_schedule_waits(bs);
    plan_waits(plan, bs);
    plan_lateness(plan);

    int late = -1;
    for (int i = plan->count - 1; i >= 0 && late < 0; --i)
        if (plan->trip[i].lateness > LATE_TOLERANCE)
            late = i;
    if (late < 0)
        return true;

    const Trip culprit = plan->trip[late];
    /* The trips of the berth on the late day and on the days that carried waiting time into it. */
    for (int n = 0, day = culprit.period; n < P; n++, day = day == 0 ? P - 1 : day - 1) {
        for (int i = 0; i < plan->count; i++)
            if (plan->trip[i].period == day && plan->trip[i].berth == culprit.berth)
                plan->trip[i].conflict = true;
        if (*bs_wait(bs, day, culprit.berth) <= LATE_TOLERANCE)
            break;
    }
    /* Every trip of the late vessel. */
    for (int i = 0; i < plan->count; i++)
        if (plan->trip[i].type == culprit.type && plan->trip[i].vessel == culprit.vessel)
            plan->trip[i].conflict = true;
    return false;
}

void plan_compact_positions(Plan *plan)
{
    for (int r = 0; r < plan->count; r++) {
        Trip *t = &plan->trip[r];
        while (t->position > 0) {
            bool taken = false;
            for (int q = 0; q < plan->count && !taken; q++)
                taken = q != r && plan->trip[q].period == t->period && plan->trip[q].berth == t->berth &&
                        plan->trip[q].position == t->position - 1;
            if (taken)
                break;
            t->position--;
        }
    }
}

/* qsort comparator: orders trips by day, berth and departing position (the timetable order). */
static int compare_trips(const void *a, const void *b)
{
    const Trip *x = a, *y = b;
    if (x->period != y->period)
        return x->period - y->period;
    if (x->berth != y->berth)
        return x->berth - y->berth;
    return x->position - y->position;
}

void plan_print(const Plan *plan, FILE *out)
{
    const Instance *in = plan->inst;
    double fixed = 0;
    for (int k = 0; k < in->num_types; k++)
        fixed += (double)plan->fleet[k] * in->type[k].fixed_cost;
    const double total = plan_cost(plan);
    fprintf(out, "objective %.6f (fleet %.6f + sailing %.6f)\n", total / COST_SCALE, fixed / COST_SCALE,
            (total - fixed) / COST_SCALE);
    fprintf(out, "fleet:");
    for (int k = 0; k < in->num_types; k++)
        fprintf(out, " type %d x %d%s", in->type[k].id, plan->fleet[k], k + 1 < in->num_types ? "," : "\n");

    Trip *sorted = XMALLOC_ARRAY(Trip, imax(plan->count, 1));
    memcpy(sorted, plan->trip, sizeof(Trip) * (size_t)plan->count);
    qsort(sorted, (size_t)plan->count, sizeof(Trip), compare_trips);
    fprintf(out, "voyages (day: berth/position vessel type#n | load m2 | wait + loading + sailing days | route):\n");
    for (int i = 0; i < plan->count; i++) {
        const Trip *t = &sorted[i];
        const Route *r = &plan->pool->route[t->route];
        const float load_days = (float)r->load_time / SECONDS_PER_DAY;
        fprintf(out, "  day %d: berth %d/%d vessel %d#%d | %7.2f | %.2f + %.2f + %.2f | berth", t->period + 1,
                t->berth >= 0 ? in->berth[t->berth].id : 0, t->position + 1, in->type[t->type].id, t->vessel + 1,
                r->load / 1e3, t->wait, load_days, (float)r->time / SECONDS_PER_DAY - load_days);
        const int *s = route_stops(plan->pool, t->route);
        for (int k = 0; k < r->size; k++)
            fprintf(out, " -> %d", in->inst[s[k]].id);
        fprintf(out, " -> berth\n");
    }
    fprintf(out, "%d voyages\n", plan->count);
    free(sorted);
}
