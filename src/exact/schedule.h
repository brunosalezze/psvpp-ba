// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file schedule.h
 * @brief Combinatorial algorithms for the scheduling problems of the
 *        feasibility subproblem (paper, Section 3): given voyages with fixed
 *        departure days and vessel types, find a vessel, a berth and a
 *        departing position for each of them.
 *
 * The paper solves these problems as MIP models. Here they are solved by
 * dynamic programming over the days of the week and by a small search, with
 * exact integer arithmetic (times in seconds); see docs/algorithms.md.
 *
 * The rules being checked (complete model, Section 2.1):
 *  - **Vessels.** A vessel is busy from the departure day of a voyage until
 *    the voyage ends, and cannot depart again before the following day. Every
 *    voyage is performed by the same vessel every week.
 *  - **Berths.** A berth loads at most ScheduleProblem::positions voyages per
 *    day, one after the other. Loading that does not fit in the day is carried
 *    over to the next day as waiting time, at most ScheduleProblem::max_carry,
 *    and the week is cyclic. Every voyage uses the same berth every week.
 *  - **Coupling.** A voyage that waits for its berth returns later: the wait
 *    plus the voyage time must end before the next departure of its vessel.
 */
#ifndef PSVPP_EXACT_SCHEDULE_H
#define PSVPP_EXACT_SCHEDULE_H

#include <stdint.h>

/** @brief Value of a waiting-time limit meaning "no limit". */
enum { WAIT_UNLIMITED = INT32_MAX };

/** @brief The voyages of a plan, reduced to what their scheduling depends on. Times are in seconds. */
typedef struct ScheduleProblem {
    int periods;      /**< days of the cyclic horizon */
    int berths;       /**< identical berths */
    int positions;    /**< departures per berth and day */
    int day_length;   /**< seconds per day */
    int max_carry;    /**< waiting time that a berth may carry into a day; less than a day */
    int num_types;    /**< vessel types */
    const int *fleet; /**< [num_types] vessels of each type */
    int n;            /**< voyages */
    const int *day;   /**< [n] departure day, 0 .. periods-1 */
    const int *type;  /**< [n] vessel type */
    const int *load;  /**< [n] loading time at the berth */
    const int *time;  /**< [n] loading plus sailing time: the vessel is busy on time / day_length + 1 days */
} ScheduleProblem;

/** @brief Answer of schedule_plan(). */
typedef enum ScheduleStatus {
    SCHEDULE_FEASIBLE,
    SCHEDULE_INFEASIBLE,
    SCHEDULE_UNDECIDED, /**< the search limit was reached before a proof (not observed on the benchmark instances) */
} ScheduleStatus;

/** @brief Reusable working memory of the algorithms (one per thread). */
typedef struct Scheduler Scheduler;

/** @brief Creates a workspace; release it with scheduler_free(). */
[[nodiscard]] Scheduler *scheduler_new(void);

/** @brief Releases a workspace (accepts @c nullptr). */
void scheduler_free(Scheduler *s);

/**
 * @brief Vessel assignment: colours circular arcs with @p fleet colours.
 *
 * Arc @c i occupies the days `start[i] .. start[i] + length[i] - 1` (cyclically)
 * of one vessel. The arcs that contain the least crowded day boundary go to
 * distinct vessels; the others are assigned in order of their first day by a
 * depth-first search that remembers the vessel states already found hopeless.
 *
 * @param s       workspace.
 * @param periods days of the cycle.
 * @param n       number of arcs.
 * @param start   [n] first day of each arc.
 * @param length  [n] days of each arc, 1 .. periods.
 * @param fleet   vessels available.
 * @param vessel  [n] receives the vessel of each arc, 0 .. fleet-1 (may be @c nullptr).
 * @return true if the arcs can be assigned.
 */
bool schedule_vessels(Scheduler *s, int periods, int n, const int *start, const int *length, int fleet, int *vessel);

/**
 * @brief Berth assignment with limited waiting times.
 *
 * Finds a berth and a departing position for the voyages such that no berth
 * carries more than @c max_carry into any day and no voyage waits longer than
 * its limit. Dynamic programme over the days: a state is the multiset, over the
 * berths, of the waiting time carried into the day; states that are no better
 * than another one on every berth are dropped. The week being cyclic, the
 * carry into the first day is not known in advance: each berth also records
 * what the days so far would allow it to be, and the last day must deliver it
 * (see the implementation notes in schedule.c).
 *
 * @param s        workspace.
 * @param p        problem (the vessel data are not used).
 * @param include  [n] voyages to schedule, or @c nullptr for all.
 * @param max_wait [n] largest waiting time of each voyage (::WAIT_UNLIMITED for none), or @c nullptr for no limit.
 * @param berth    [n] receives the berth of each scheduled voyage (may be @c nullptr).
 * @param position [n] receives its departing position, 0 = first of the day (may be @c nullptr).
 * @return true if a schedule exists.
 */
bool schedule_berths(Scheduler *s, const ScheduleProblem *p, const bool *include, const int *max_wait, int *berth,
                     int *position);

/**
 * @brief Vessels and berths together: the feasibility subproblem.
 *
 * A voyage that waits @c w for its berth blocks its vessel for
 * `ceil((w + time) / day_length)` days, i.e. 0, 1 or 2 days more than without
 * waiting. The search looks for a number of extra days per voyage for which the
 * vessels can still be assigned (schedule_vessels() with longer arcs) and the
 * berths can keep every wait within the corresponding limit
 * (schedule_berths()). It starts without extra days; when the berths fail, it
 * identifies a minimal set of voyages whose limits cause the failure and tries
 * to give one more day to each of them in turn.
 *
 * @param s        workspace.
 * @param p        problem.
 * @param vessel   [n] receives the vessel of each voyage within its type (may be @c nullptr).
 * @param berth    [n] receives its berth (may be @c nullptr).
 * @param position [n] receives its departing position (may be @c nullptr).
 */
ScheduleStatus schedule_plan(Scheduler *s, const ScheduleProblem *p, int *vessel, int *berth, int *position);

/**
 * @brief Checks a complete assignment directly from the rules (independent of the algorithms above).
 * @return true if the assignment is feasible.
 */
bool schedule_check(const ScheduleProblem *p, const int *vessel, const int *berth, const int *position);

/** @brief Work counters of a workspace since its creation (for statistics and tests). */
typedef struct ScheduleStats {
    long berth_calls;  /**< schedule_berths() calls */
    long berth_states; /**< day-to-day transitions examined by the berth search */
    long berth_limits; /**< schedule_berths() calls stopped at the work limit */
    long vessel_calls; /**< schedule_vessels() calls */
    long vessel_nodes; /**< nodes of the vessel search */
    long plan_calls;   /**< schedule_plan() calls */
    long plan_nodes;   /**< nodes of the search of schedule_plan() */
    long plan_waits;   /**< feasible plans in which some vessel is blocked longer because of a wait */
    long plan_limits;  /**< schedule_plan() calls stopped at a limit (answer ::SCHEDULE_UNDECIDED) */
} ScheduleStats;

/** @brief Returns the counters of @p s. */
ScheduleStats scheduler_stats(const Scheduler *s);

/**
 * @brief true if the last schedule_berths() call stopped at its work limit: its answer `false`
 *        then means "not found", not "none exists".
 */
bool scheduler_gave_up(const Scheduler *s);

/**
 * @brief After schedule_berths() returned false (without giving up): the voyages departing on
 *        the @p days consecutive days from day @p first (cyclically) cannot be scheduled on their
 *        own, with the limits that were given.
 */
void scheduler_conflict_days(const Scheduler *s, int *first, int *days);

#endif /* PSVPP_EXACT_SCHEDULE_H */
