// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file moves.h
 * @brief Incremental evaluation and application of segment insertions and removals.
 *
 * Every modification of a heuristic solution is expressed as inserting or
 * removing a sequence of consecutive visits (a @e segment) between two nodes of
 * one voyage. The evaluator computes the exact change of the penalised
 * objective in time independent of the solution size by propagating the
 * consequences of the move:
 *
 *  1. the voyage's distance, load, loading time and sailing time;
 *  2. the fleet: opening a new voyage may activate a vessel (fixed cost) and
 *     splits the time window of the vessel's previous voyage; closing the last
 *     voyage releases it;
 *  3. the berth: the voyage's loading time changes the occupation of its berth
 *     on its departure day; occupation beyond one day is carried over to the
 *     next days as waiting time (cyclically), and occupation beyond the
 *     tolerance TOL (or beyond the week) is a berth time violation;
 *  4. the waiting times of the voyages queued after it at the berth, which may
 *     make their vessels late for their next departures.
 */
#ifndef PSVPP_HEURISTIC_MOVES_H
#define PSVPP_HEURISTIC_MOVES_H

#include "heuristic/search.h"
#include "heuristic/solution.h"

#include <math.h>

/** @brief Direction of a segment move (the sign of every incremental quantity). */
typedef enum MoveDir { MOVE_REMOVE = -1, MOVE_INSERT = +1 } MoveDir;

/**
 * @brief Evaluates, and optionally applies, a segment move.
 *
 * @param S       search context (penalty weights and scratch buffers).
 * @param s       solution.
 * @param seg     the visit nodes of the segment, in route order.
 * @param n       segment length (>= 1).
 * @param before  node preceding the segment (a voyage head or a visit).
 * @param after   node following the segment; for an insertion it must be
 *                `s->node[before].next`, for a removal the node after the
 *                last segment node.
 * @param dir     ::MOVE_INSERT to insert the (unrouted) segment between
 *                @p before and @p after, ::MOVE_REMOVE to remove it.
 * @param commit  when false only the cost change is computed; when true the
 *                solution is updated accordingly (the penalised cost field
 *                Solution::cost is @e not updated; see sol_update_cost()).
 * @return Change of the penalised objective caused by the move.
 */
Cost move_segment(Search *S, Solution *s, const int *seg, int n, int before, int after, MoveDir dir, bool commit);

/** @brief Cost of removing visit @p node from its voyage (no change is made). */
static inline Cost move_remove_cost(Search *S, Solution *s, int node)
{
    return move_segment(S, s, &node, 1, s->node[node].prev, s->node[node].next, MOVE_REMOVE, false);
}

/** @brief Cost of inserting the unrouted visit @p node right after @p before (no change is made). */
static inline Cost move_insert_cost(Search *S, Solution *s, int node, int before)
{
    return move_segment(S, s, &node, 1, before, s->node[before].next, MOVE_INSERT, false);
}

/**
 * @brief Change of the un-penalised cost of a vessel of type @p t whose sailed
 *        distance changes by @p dist and that is activated (@p vessel_event =
 *        ::MOVE_INSERT), released (::MOVE_REMOVE) or neither (0).
 */
static inline int vessel_cost_delta(const VesselType *t, int dist, int vessel_event)
{
    return (int)((float)dist * t->var_cost + (float)(vessel_event * t->fixed_cost));
}

/**
 * @brief Second part of move_insert_cost_below(): evaluates the insertion of
 *        @p node after @p before, whose detour is @p detour, stopping as soon
 *        as the cost is known to reach @p bound (INFINITY = never stop).
 */
Cost move_insert_cost_bounded(Search *S, Solution *s, int node, int before, int detour, Cost bound);

/**
 * @brief Cost of inserting the unrouted visit @p node right after @p before, when it is below @p bound.
 *
 * Insertions are evaluated to find the cheapest ones, so the exact cost of a
 * position that cannot compete is not needed. An insertion never reduces a
 * violation: it only adds load, loading time and (unless the rounding of the
 * distances makes the detour negative) sailing time, so every term of its
 * cost is non-negative and the terms computed so far are a lower bound on it.
 * The evaluation stops as soon as that lower bound reaches @p bound. The first
 * term, the cost of the detour and of the vessel that the insertion activates,
 * needs no berth or time-window computation and already discards about half of
 * the positions.
 *
 * @param bound value from which the exact cost is of no interest to the caller.
 * @return The exact cost change if it is below @p bound; otherwise a value in
 *         [@p bound, exact cost]. No change is made.
 */
static inline Cost move_insert_cost_below(Search *S, Solution *s, int node, int before, Cost bound)
{
    const Instance *in = S->inst;
    const VoyageKey at = s->node[before].at;
    const int after = s->node[before].next;
    const int detour =
        inst_node_dist(in, before, node) + inst_node_dist(in, node, after) - inst_node_dist(in, before, after);
    if (detour < 0)
        return move_insert_cost_bounded(S, s, node, before, detour, INFINITY);
    const int activates = s->vessel[sol_vessel_of(s, at)].visits == 0 ? MOVE_INSERT : 0;
    const Cost floor = (float)vessel_cost_delta(&in->type[at.type], detour, activates);
    return floor >= bound ? floor : move_insert_cost_bounded(S, s, node, before, detour, bound);
}

/** @brief Removes visit @p node from its voyage and returns the cost change. */
static inline Cost move_remove(Search *S, Solution *s, int node)
{
    return move_segment(S, s, &node, 1, s->node[node].prev, s->node[node].next, MOVE_REMOVE, true);
}

/** @brief Inserts the unrouted visit @p node right after @p before and returns the cost change. */
static inline Cost move_insert(Search *S, Solution *s, int node, int before)
{
    return move_segment(S, s, &node, 1, before, s->node[before].next, MOVE_INSERT, true);
}

/** @brief Inserts @p n unrouted visits, in order, right after @p before. */
static inline Cost move_insert_many(Search *S, Solution *s, const int *nodes, int n, int before)
{
    return move_segment(S, s, nodes, n, before, s->node[before].next, MOVE_INSERT, true);
}

#endif /* PSVPP_HEURISTIC_MOVES_H */
