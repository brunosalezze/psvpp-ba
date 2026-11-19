// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file plan.h
 * @brief A set of scheduled voyages (an integer solution of the exact models)
 *        and the berth/vessel schedule computations used to check it.
 *
 * Times in this module are expressed in days (periods), as in the models.
 */
#ifndef PSVPP_EXACT_PLAN_H
#define PSVPP_EXACT_PLAN_H

#include "common/instance.h"
#include "exact/routes.h"

#include <stdint.h>
#include <stdio.h>

/** @brief Tolerance TOL of the berth occupation per day, in days (1.5). */
static inline double tol_days(void) { return (double)BERTH_TOLERANCE / SECONDS_PER_DAY; }

/** @brief Waiting time that a day may carry into the next one, TOL - 1, in days (0.5). */
static inline double max_carry_days(void) { return tol_days() - 1.0; }

/** @brief One voyage of a plan with its (partial) assignment. */
typedef struct Trip {
    long route;      /**< voyage in the RoutePool */
    int8_t period;   /**< departure day t */
    int8_t type;     /**< vessel type v */
    int8_t vessel;   /**< vessel n of type v (-1 = not assigned) */
    int8_t berth;    /**< berth b (-1 = not assigned) */
    int8_t position; /**< departing position k at the berth (-1 = not assigned) */
    bool conflict;   /**< member of the conflict set being built */
    float fat_time;  /**< slack until the vessel's next departure (days) */
    float wait;      /**< waiting time for the berth (days) */
    float big_m;     /**< big-M of the waiting-time constraint of the subproblem */
    float lateness;  /**< max(wait - fat_time, 0) */
    int var;         /**< MIP column that selected the trip */
} Trip;

/** @brief A set of trips with the fleet that performs them. */
typedef struct Plan {
    const Instance *inst;
    const RoutePool *pool;
    Trip *trip;   /**< [count] */
    int count;    /**< trips in use */
    int capacity; /**< allocated trips */
    int *fleet;   /**< [num_types] vessels of each type */
} Plan;

/** @brief Waiting and loading times of every berth and day (days). */
typedef struct BerthSchedule {
    const Instance *inst;
    float *wait; /**< [num_periods * num_berths] time carried over from the previous day, \f$t^{WB}\f$ */
    float *load; /**< [num_periods * num_berths * MAX_DEPARTURES_PER_BERTH] loading time per position, \f$t^{LB}\f$ */
} BerthSchedule;

/** @name Plans */
///@{
void plan_init(Plan *plan, const Instance *inst, const RoutePool *pool);
void plan_free(Plan *plan);
void plan_clear(Plan *plan);
void plan_copy(Plan *dst, const Plan *src);
void plan_add(Plan *plan, Trip trip);
/** @brief Fixed plus sailing cost of the plan (internal cost units). */
double plan_cost(const Plan *plan);
/** @brief Loading time of the voyage of a trip (days). */
float trip_load_days(const Plan *plan, const Trip *t);
///@}

/** @name Schedules */
///@{
void berth_schedule_init(BerthSchedule *bs, const Instance *inst);
void berth_schedule_free(BerthSchedule *bs);
static inline float *bs_wait(const BerthSchedule *bs, int period, int berth)
{
    return &bs->wait[period * bs->inst->num_berths + berth];
}
static inline float *bs_load(const BerthSchedule *bs, int period, int berth, int pos)
{
    return &bs->load[(period * bs->inst->num_berths + berth) * MAX_DEPARTURES_PER_BERTH + pos];
}

/**
 * @brief Waiting times implied by the loading times: the occupation of a
 *        berth beyond one day is carried over to the next day, cyclically
 *        (constraints 7-8 of the complete model taken as equalities).
 */
void berth_schedule_waits(BerthSchedule *bs);

/**
 * @brief Slack of every trip: days until the next departure of the same
 *        vessel minus the voyage duration (requires the vessel assignment).
 */
void plan_fat_times(Plan *plan);

/** @brief Waiting time of every trip given the berth schedule (requires berth and position). */
void plan_waits(Plan *plan, const BerthSchedule *bs);

/** @brief Lateness of every trip: max(wait - fat_time, 0). */
void plan_lateness(Plan *plan);

/**
 * @brief Big-M of the waiting-time constraint of every trip in the subproblem model: TOL - 1 plus
 *        the longest loading time of its day.
 *
 * The trip itself counts: the constraint of the second position also binds a trip loading in
 * the first position, with its own loading time on the left-hand side. (The original
 * implementation took the longest loading time of the *other* trips, which could declare a
 * feasible plan infeasible.)
 */
void plan_big_m(Plan *plan);

/**
 * @brief Checks a fully assigned plan for vessels that cannot leave on time.
 *
 * When some trip is late, the conflict flags mark the trips responsible: the
 * trips of the late trip's berth on its day and on the preceding days that
 * carried waiting time into it, and every trip of the late vessel.
 *
 * @return true if no trip is late.
 */
bool plan_check_overlaps(Plan *plan, BerthSchedule *bs);
///@}

/**
 * @brief Moves voyages to earlier departing positions left empty at their
 *        berth and day (the loading order, hence the schedule, is unchanged).
 */
void plan_compact_positions(Plan *plan);

/** @brief Prints the plan (fleet, cost and the timetable of every trip). */
void plan_print(const Plan *plan, FILE *out);

#endif /* PSVPP_EXACT_PLAN_H */
