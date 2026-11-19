// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file moves.c
 * @brief Incremental segment insertion/removal (see moves.h).
 *
 * The evaluation is split into the phases listed in moves.h. All the deltas
 * are accumulated in a ::MoveDelta; commit_move() then applies them. Integer
 * truncations and the float evaluation order are those of the original code,
 * so that removals exactly undo insertions.
 */
#include "heuristic/moves.h"

#include "common/util.h"

#include <math.h>
#include <string.h>

/** @brief All the consequences of one segment move. */
typedef struct MoveDelta {
    /* where */
    VoyageKey at;     /**< voyage receiving/losing the segment */
    int vessel;       /**< its vessel slot */
    int queue_pos;    /**< its position in the departure queue (the end of the queue for a new voyage) */
    int voyage_event; /**< MOVE_INSERT: opens the voyage; MOVE_REMOVE: empties it; 0 otherwise */
    int vessel_event; /**< MOVE_INSERT: activates the vessel; MOVE_REMOVE: deactivates it; 0 otherwise */
    /* changes of the voyage (dist is set by the caller, the others by the phases) */
    int dist, load, load_time, travel, wait, total;
    int visits;   /**< change of the number of visits of the voyage (signed) */
    int max_time; /**< time window of the voyage (taken from the previous voyage when opening one) */
    /* changes of costs and violations */
    int cost;         /**< fixed + sailing cost */
    int excess_load;  /**< capacity violation of the voyage */
    int excess_time;  /**< time violation of the vessel (this voyage and its predecessor) */
    int berth_excess; /**< berth occupation beyond TOL per day and beyond the week */
    int queue_excess; /**< time violations of the voyages queued after this one */
    int d_departures; /**< departures beyond the berth limit */
    int num_queue;    /**< entries of Search::queue in use */
    bool carries;     /**< the move changes the time carried between days (Search::shift is not all zero) */
} MoveDelta;

/* Distance sailed with the segment between @p before and @p after minus the direct leg between them. */
static int segment_detour(const Instance *in, const int *seg, int n, int before, int after)
{
    int dist =
        inst_node_dist(in, before, seg[0]) + inst_node_dist(in, seg[n - 1], after) - inst_node_dist(in, before, after);
    for (int i = 0; i + 1 < n; i++)
        dist += inst_node_dist(in, seg[i], seg[i + 1]);
    return dist;
}

/* Capacity violation added by the move: load beyond the vessel capacity, plus Search::size_excess for
 * every installation beyond the largest voyage size. Both are weighed by PenaltyWeights::load. */
static int load_excess_delta(const Search *S, const Voyage *r, const VesselType *tc, const MoveDelta *d)
{
    return excess_delta(r->load - tc->capacity, d->load) +
           S->size_excess * excess_delta(r->size - S->par->max_route_size, d->visits);
}

/* Phase 1: load, times of the segment and the voyage-opening wait (MoveDelta::dist is the signed detour). */
static void delta_voyage(const Search *S, const Solution *s, MoveDelta *d, const int *seg, int n, MoveDir dir)
{
    const Instance *in = S->inst;
    const int b = d->at.berth, p = d->at.period;
    const Voyage *r = sol_voyage(s, d->vessel, p);
    const Vessel *ve = &s->vessel[d->vessel];

    if (dir == MOVE_INSERT && r->size == 0) {
        d->voyage_event = MOVE_INSERT;
        if (ve->visits == 0)
            d->vessel_event = MOVE_INSERT;
    }
    if (dir == MOVE_REMOVE && r->size == n) {
        d->voyage_event = MOVE_REMOVE;
        if (ve->visits == n)
            d->vessel_event = MOVE_REMOVE;
    }

    const int dist = dir * d->dist; /* unsigned detour */
    int travel = 0, load = 0, load_time = 0;
    for (int i = 0; i < n; i++) {
        travel += in->visit_time[seg[i]];
        load += in->visit_load[seg[i]];
        load_time += in->visit_load_time[seg[i] * in->num_berths + b];
    }
    travel = (int)((float)travel + in->type[d->at.type].time_per_dist * (float)dist);

    int wait = 0;
    if (d->voyage_event != 0) {
        /* Opening (closing) a voyage adds (removes) the berth set-up time and the
         * waiting time before its turn in the queue. */
        const BerthDay *day = sol_day(s, b, p);
        load_time += in->berth[b].setup_time;
        wait = day->carry_time;
        if (d->voyage_event == MOVE_INSERT) {
            wait += day->load_time;
        } else {
            const VoyageKey *q = sol_queue(s, b, p);
            for (int i = d->queue_pos - 1; i >= 0; --i)
                wait += sol_voyage(s, inst_vessel(in, b, q[i].type, q[i].slot), p)->load_time;
        }
    }
    d->load = dir * load;
    d->load_time = dir * load_time;
    d->travel = dir * travel;
    d->wait = dir * wait;
}

/* Phase 2: time window of the voyage and the departure limit of the berth. */
static void delta_time_window(const Search *S, const Solution *s, MoveDelta *d, MoveDir dir)
{
    const Instance *in = S->inst;
    const int p = d->at.period;
    const Voyage *r = sol_voyage(s, d->vessel, p);
    if (dir == MOVE_INSERT && r->size == 0) {
        if (sol_day(s, d->at.berth, p)->num_departures >= MAX_DEPARTURES_PER_BERTH)
            d->d_departures++;
        if (s->vessel[d->vessel].visits == 0) {
            d->max_time = in->horizon;
        } else {
            /* The new voyage takes the end of the previous voyage's window. */
            int gap = p - r->prev_period;
            if (gap < 0)
                gap += in->num_periods;
            d->max_time = sol_voyage(s, d->vessel, r->prev_period)->max_time - gap * SECONDS_PER_DAY;
        }
    } else {
        d->max_time = r->max_time;
    }
}

/* Phase 3: propagation of the berth occupation over the following days.
 * A berth occupied for more than one day carries the excess to the next day
 * (cyclically). When the weekly loading time exceeds the week, the overflow is
 * absorbed by the 'sink' of the days so that the propagation terminates. */
static void delta_berth(Search *S, const Solution *s, MoveDelta *d)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, p = d->at.period;
    const BerthLog *bl = &s->berth[d->at.berth];
    const BerthDay *days = sol_day(s, d->at.berth, 0);
    int *shift = S->shift, *dsink = S->d_sink;
    if (S->shift_dirty) { /* both buffers are all zero otherwise */
        memset(shift, 0, sizeof(int) * (size_t)P);
        memset(dsink, 0, sizeof(int) * (size_t)P);
        S->shift_dirty = false;
    }

    int overflow = excess_delta(bl->total_load_time - in->horizon, d->load_time);
    int dsk = ipos(days[p].sink + overflow) - days[p].sink;
    int used = days[p].carry_time + days[p].load_time - days[p].sink;
    int new_used = used + d->load_time - dsk;
    int carry = ipos(new_used - SECONDS_PER_DAY) - ipos(used - SECONDS_PER_DAY);
    d->berth_excess = ipos(new_used - BERTH_TOLERANCE) - ipos(used - BERTH_TOLERANCE) + overflow;
    overflow -= dsk;

    d->carries = carry != 0 || overflow != 0;
    if (dsk != 0 || d->carries) {
        dsink[p] = dsk;
        S->shift_dirty = true;
    }
    for (int i = p, steps = 0; carry != 0 || overflow != 0;) {
        if (++i >= P)
            i = 0;
        const BerthDay *day = &days[i];
        dsk = ipos(day->sink + dsink[i] + overflow) - day->sink - dsink[i];
        used = day->carry_time + day->load_time - day->sink - dsink[i] + shift[i];
        if (i == p)
            used += d->load_time;
        new_used = used + carry - dsk;
        overflow -= dsk;
        dsink[i] += dsk;
        shift[i] += carry;
        carry = ipos(new_used - SECONDS_PER_DAY) - ipos(used - SECONDS_PER_DAY);
        d->berth_excess += ipos(new_used - BERTH_TOLERANCE) - ipos(used - BERTH_TOLERANCE);
        if (++steps > 2 * P)
            break;
    }
}

/* Records that the queued voyage @p k of day @p period waits @p diff seconds longer. */
static void queue_shift(Search *S, const Solution *s, MoveDelta *d, int period, VoyageKey k, int diff)
{
    const int slack = voyage_time_slack(sol_voyage(s, inst_vessel(S->inst, d->at.berth, k.type, k.slot), period));
    const int excess = excess_delta(slack, diff);
    d->queue_excess += excess;
    S->queue[d->num_queue++] =
        (QueueShift){.period = (int8_t)period, .type = k.type, .slot = k.slot, .d_wait = diff, .d_excess = excess};
}

/* Phase 4: waiting-time changes of the voyages queued at the berth. */
static void delta_queue(Search *S, const Solution *s, MoveDelta *d)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, p = d->at.period, b = d->at.berth;
    const Voyage *r = sol_voyage(s, d->vessel, p);
    d->queue_excess = 0;
    d->num_queue = 0;
    if (!d->carries) {
        /* Usual case: no other day is affected, so only the voyages queued behind this one on its
         * own day wait for its loading-time change. */
        if (d->load_time != 0) {
            const int num = sol_day(s, b, p)->num_departures;
            const VoyageKey *q = sol_queue(s, b, p);
            for (int j = d->queue_pos + 1; j < num; j++)
                queue_shift(S, s, d, p, q[j], d->load_time);
        }
        return;
    }
    if (d->voyage_event != MOVE_REMOVE)
        d->wait += S->shift[p];
    for (int i = P - 1; i >= 0; --i) {
        if (S->shift[i] == 0 && i != p)
            continue;
        int behind = 0; /* loading-time change of the voyages ahead in the queue */
        const int num = sol_day(s, b, i)->num_departures;
        const VoyageKey *q = sol_queue(s, b, i);
        for (int j = 0; j < num; j++) {
            if (i == p && j == d->queue_pos) {
                behind += d->load_time;
                continue;
            }
            if (S->shift[i] + behind != 0)
                queue_shift(S, s, d, i, q[j], S->shift[i] + behind);
        }
    }
    /* The previous voyage of the same vessel is accounted for in delta_costs(). */
    if (d->voyage_event != 0 && r->prev_period >= 0 && r->prev_period != p)
        for (int i = d->num_queue - 1; i >= 0; --i) {
            QueueShift *e = &S->queue[i];
            if (e->type == d->at.type && e->slot == d->at.slot && e->period == r->prev_period) {
                d->queue_excess -= e->d_excess;
                e->d_excess = 0;
            }
        }
}

/* Phase 5: total duration, costs and violations of the voyage. */
static void delta_costs(const Search *S, const Solution *s, MoveDelta *d)
{
    const Instance *in = S->inst;
    const VesselType *tc = &in->type[d->at.type];
    const int p = d->at.period;
    const Voyage *r = sol_voyage(s, d->vessel, p);

    d->total = d->load_time + d->travel + d->wait;
    if (d->voyage_event == MOVE_REMOVE && sol_day(s, d->at.berth, p)->num_departures > MAX_DEPARTURES_PER_BERTH)
        d->d_departures--;

    d->cost = vessel_cost_delta(tc, d->dist, d->vessel_event);
    d->excess_load = load_excess_delta(S, r, tc, d);
    d->excess_time = excess_delta(r->time - d->max_time, d->total);
    if (d->voyage_event != 0 && r->prev_period >= 0 && r->prev_period != p) {
        /* Opening (closing) a voyage shrinks (extends) the previous voyage's window. */
        const int slack = voyage_time_slack(sol_voyage(s, d->vessel, r->prev_period));
        const int shift = S->shift[r->prev_period];
        if (d->voyage_event == MOVE_INSERT)
            d->excess_time += ipos(slack + shift + d->max_time) - ipos(slack);
        else
            d->excess_time += ipos(slack + shift - d->max_time) - ipos(slack);
    }
}

/* Unlinks/links the segment and updates the voyage structure of the vessel. */
static void commit_structure(const Search *S, Solution *s, const MoveDelta *d, const int *seg, int n, int before,
                             int after, MoveDir dir)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, p = d->at.period, b = d->at.berth, v = d->vessel;
    const VoyageKey rd = d->at;
    Node *nd = s->node;
    Voyage *r = sol_voyage(s, v, p);
    Vessel *ve = &s->vessel[v];
    BerthDay *day = sol_day(s, b, p);
    BerthLog *bl = &s->berth[b];
    VoyageKey *q = sol_queue(s, b, p);
    int8_t *active = sol_active(s, b, rd.type);

    if (dir == MOVE_REMOVE) {
        nd[before].next = after;
        nd[after].prev = before;
        if (r->size == n) { /* the voyage becomes empty: leave the queue */
            if (day->num_departures-- > MAX_DEPARTURES_PER_BERTH)
                bl->excess_departures--;
            for (int i = r->queue_pos; i < day->num_departures; i++) {
                const VoyageKey k = q[i + 1];
                sol_voyage(s, inst_vessel(in, b, k.type, k.slot), p)->queue_pos--;
                q[i] = k;
            }
            r->queue_pos = -1;
            if (ve->visits == n) { /* the vessel becomes inactive */
                ve->first_period = -1;
                for (int i = 0; i < P; i++)
                    sol_voyage(s, v, i)->prev_period = -1;
                /* The active prefix ends at the last slot still in use (it may contain idle slots
                 * until sol_compact()); being a function of the slots in use, it is restored when
                 * the move is undone. */
                int last = *active - 1;
                while (last >= 0 && (last == rd.slot || sol_vessel(s, b, rd.type, last)->visits == 0))
                    last--;
                *active = (int8_t)(last + 1);
            } else { /* merge the window into the previous voyage */
                const int prv = r->prev_period, nxt = r->next_period;
                sol_voyage(s, v, prv)->max_time += r->max_time;
                sol_voyage(s, v, nxt)->prev_period = (int8_t)prv;
                sol_voyage(s, v, prv)->next_period = (int8_t)nxt;
                if (ve->first_period == p)
                    ve->first_period = nxt;
                for (int i = p;; i++) {
                    if (i >= P)
                        i -= P;
                    if (i == nxt)
                        break;
                    sol_voyage(s, v, i)->prev_period = (int8_t)prv;
                }
            }
            r->max_time = 0;
            r->next_period = -1; /* as in an untouched empty voyage: undoing an insertion is then exact */
        }
        return;
    }

    /* insertion: link seg[0..n-1] between before and after */
    nd[before].next = seg[0];
    nd[after].prev = seg[n - 1];
    nd[seg[0]].prev = before;
    nd[seg[n - 1]].next = after;
    nd[seg[0]].at = rd;
    if (n > 1) {
        nd[seg[0]].next = seg[1];
        nd[seg[n - 1]].prev = seg[n - 2];
        nd[seg[n - 1]].at = rd;
        for (int i = 1; i < n - 1; i++) {
            nd[seg[i]].prev = seg[i - 1];
            nd[seg[i]].next = seg[i + 1];
            nd[seg[i]].at = rd;
        }
    }
    if (r->size == 0) { /* a new voyage: join the queue and split the vessel's week */
        r->queue_pos = day->num_departures;
        q[day->num_departures] = rd;
        if (++day->num_departures > MAX_DEPARTURES_PER_BERTH)
            bl->excess_departures++;
        if (ve->visits == 0) {
            ve->first_period = p;
            r->next_period = r->prev_period = (int8_t)p;
            r->max_time = in->horizon;
            if (*active <= rd.slot)
                *active = (int8_t)(rd.slot + 1);
        } else {
            const int prv = r->prev_period, nxt = sol_voyage(s, v, prv)->next_period;
            r->next_period = (int8_t)nxt;
            sol_voyage(s, v, nxt)->prev_period = (int8_t)p;
            sol_voyage(s, v, prv)->next_period = (int8_t)p;
            r->max_time = d->max_time;
            sol_voyage(s, v, prv)->max_time -= d->max_time;
        }
        const int nxt = r->next_period;
        for (int i = p + 1;; i++) {
            if (i >= P)
                i -= P;
            if (i == nxt)
                break;
            sol_voyage(s, v, i)->prev_period = (int8_t)p;
        }
    }
}

/* Applies every delta of the move to the solution. */
static void commit_move(const Search *S, Solution *s, const MoveDelta *d, const int *seg, int n, int before, int after,
                        MoveDir dir)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, p = d->at.period, b = d->at.berth, N = in->num_installations;
    Voyage *r = sol_voyage(s, d->vessel, p);
    Vessel *ve = &s->vessel[d->vessel];
    BerthDay *days = sol_day(s, b, 0);
    BerthLog *bl = &s->berth[b];

    commit_structure(S, s, d, seg, n, before, after, dir);

    ve->visits += n * dir;
    r->size += n * dir;
    for (int i = 0; i < n; i++) {
        const int u = in->visit_inst[seg[i]];
        s->visit_count[p * N + u] = (int8_t)(s->visit_count[p * N + u] + dir);
        s->visit_node[p * N + u] = seg[i];
    }
    for (int i = d->num_queue - 1; i >= 0; --i) {
        const QueueShift *e = &S->queue[i];
        const int w = inst_vessel(in, b, e->type, e->slot);
        Voyage *rq = sol_voyage(s, w, e->period);
        rq->wait_time += e->d_wait;
        rq->time += e->d_wait;
        s->vessel[w].excess_time += e->d_excess;
    }
    for (int i = P - 1; i >= 0; --i) {
        days[i].carry_time += S->shift[i];
        days[i].sink += S->d_sink[i];
    }
    r->dist += d->dist;
    r->load += d->load;
    r->load_time += d->load_time;
    r->wait_time += d->wait;
    r->travel_time += d->travel;
    r->time += d->total;
    days[p].load_time += d->load_time;
    bl->total_load_time += d->load_time;
    bl->excess_time += d->berth_excess;
    ve->excess_load += d->excess_load;
    ve->excess_time += d->excess_time;
    ve->cost += d->cost;
}

/* Runs the phases for the move described by @p d (voyage, vessel and signed detour already set).
 *
 * @p bound is INFINITY except for an insertion with a non-negative detour whose cost is only needed
 * when it is below @p bound. For such a move every term of the cost is non-negative and the berth
 * (phases 3-4) can only add waiting time, so the terms known after phase 2, and then after phase 3,
 * are lower bounds on the cost; when one of them reaches @p bound it is returned instead of the cost. */
static Cost evaluate(Search *S, Solution *s, MoveDelta *d, const int *seg, int n, int before, int after, MoveDir dir,
                     Cost bound, bool commit)
{
    const Voyage *r = sol_voyage(s, d->vessel, d->at.period);
    const PenaltyWeights w = S->w;
    d->queue_pos =
        (dir == MOVE_INSERT && r->size == 0) ? sol_day(s, d->at.berth, d->at.period)->num_departures : r->queue_pos;

    d->visits = dir * n;
    delta_voyage(S, s, d, seg, n, dir);
    delta_time_window(S, s, d, dir);
    if (bound < INFINITY) {
        const VesselType *tc = &S->inst->type[d->at.type];
        const int time_excess = excess_delta(r->time - d->max_time, d->load_time + d->travel + d->wait);
        const int load_excess = load_excess_delta(S, r, tc, d);
        const Cost partial = (float)vessel_cost_delta(tc, d->dist, d->vessel_event) + (float)time_excess * w.time +
                             (float)load_excess * w.load + (float)d->d_departures * w.berth;
        if (partial >= bound)
            return partial;
        delta_berth(S, s, d);
        if (partial + (float)d->berth_excess * w.time >= bound)
            return partial + (float)d->berth_excess * w.time;
    } else {
        delta_berth(S, s, d);
    }
    delta_queue(S, s, d);
    delta_costs(S, s, d);

    const Cost value = (float)d->cost + (float)d->excess_time * w.time + (float)d->excess_load * w.load +
                       (float)d->d_departures * w.berth + (float)d->berth_excess * w.time +
                       (float)d->queue_excess * w.time;
    if (commit)
        commit_move(S, s, d, seg, n, before, after, dir);
    return value;
}

Cost move_segment(Search *S, Solution *s, const int *seg, int n, int before, int after, MoveDir dir, bool commit)
{
    MoveDelta d = {.at = s->node[before].at};
    d.vessel = sol_vessel_of(s, d.at);
    d.dist = dir * segment_detour(S->inst, seg, n, before, after);
    return evaluate(S, s, &d, seg, n, before, after, dir, INFINITY, commit);
}

Cost move_insert_cost_bounded(Search *S, Solution *s, int node, int before, int detour, Cost bound)
{
    MoveDelta d = {.at = s->node[before].at, .dist = detour};
    d.vessel = sol_vessel_of(s, d.at);
    return evaluate(S, s, &d, &node, 1, before, s->node[before].next, MOVE_INSERT, bound, false);
}
