// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file solution.h
 * @brief Solution representation of the ALNS heuristic (paper, Section 4).
 *
 * ### Structure
 * The fleet is a fixed set of @e vessel slots (Instance::num_vessels): each
 * berth hosts Instance::type[k].max_slots slots of every type @c k. A slot is
 * @e active when it performs at least one voyage; per (berth, type) the slots
 * `0 .. active-1` up to the last active one form the active prefix, which
 * sol_compact() keeps free of idle slots between iterations.
 *
 * Every (slot, day) pair owns one ::Voyage, stored as a circular doubly-linked
 * list of nodes that starts and ends at the voyage's head node
 * (inst_head_node()). An empty voyage is a head node linked to itself.
 *
 * For every berth and day a ::BerthDay records the departures queued at the
 * berth (their order is the loading order), the loading time and the waiting
 * time carried over from the previous day, which together give the berth
 * schedule of the paper's complete model.
 *
 * All the arrays of a solution live in one memory block whose layout depends
 * only on the instance, so sol_copy() is a single @c memcpy.
 *
 * ### Costs and feasibility
 * All quantities are integers in the internal units of instance.h. Constraint
 * violations are kept as @e excess amounts (time or load beyond the limit) that
 * the penalised objective weighs with the adaptive penalty weights of the
 * ::Search (paper, Section 4.2).
 */
#ifndef PSVPP_HEURISTIC_SOLUTION_H
#define PSVPP_HEURISTIC_SOLUTION_H

#include "common/instance.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/** @brief Floating-point type of (penalised) objective values. */
typedef float Cost;

/** @brief Link of a node in the circular list of its voyage. */
typedef struct Node {
    int next;     /**< following node of the voyage (the head after the last visit) */
    int prev;     /**< preceding node of the voyage */
    VoyageKey at; /**< voyage the node currently belongs to (meaningful for visits in a voyage) */
} Node;

/** @brief Voyage of one vessel slot departing on one day. */
typedef struct Voyage {
    int size;           /**< number of visits */
    int queue_pos;      /**< position in the berth's departure queue of its day (-1 when empty) */
    int8_t next_period; /**< next day on which the same vessel departs (circular; -1 when empty) */
    int8_t prev_period; /**< previous departure day of the vessel (kept for every day, used or not) */
    int dist;           /**< sailed distance */
    int load;           /**< cargo loaded */
    int time;           /**< total duration: wait + load_time + travel_time (s) */
    int max_time;       /**< time available until the vessel's next departure (s) */
    int travel_time;    /**< sailing and service time at the installations (s) */
    int wait_time;      /**< waiting time for the berth (s) */
    int load_time;      /**< loading time at the berth, including the set-up time (s) */
} Voyage;

/** @brief Aggregated data of one vessel slot over the week. */
typedef struct Vessel {
    int visits;       /**< number of visits over all its voyages (0 = inactive slot) */
    int first_period; /**< first departure day (-1 when inactive) */
    int excess_time;  /**< total time violation of its voyages (overlaps with the next departure) */
    int excess_load;  /**< total capacity violation of its voyages */
    int cost;         /**< un-penalised cost: fixed cost (if active) plus sailing costs */
} Vessel;

/** @brief State of one berth on one day. */
typedef struct BerthDay {
    int num_departures; /**< voyages departing on this day (queue length) */
    int load_time;      /**< total loading time of the departures (s) */
    int carry_time;     /**< waiting time imposed by the previous day, \f$t^{WB}\f$ (s) */
    int sink;           /**< share of the carried time absorbed by the weekly overflow (s) */
} BerthDay;

/** @brief Weekly aggregates of one berth. */
typedef struct BerthLog {
    int total_load_time;   /**< loading time over the week (s) */
    int excess_time;       /**< berth time violation (beyond TOL per day and beyond the week) */
    int excess_departures; /**< departures beyond ::MAX_DEPARTURES_PER_BERTH, summed over the days */
} BerthLog;

/** @brief Number of each kind of infeasibility found by sol_check_feasibility(). */
typedef struct Violations {
    int load;  /**< vessels exceeding their capacity */
    int time;  /**< vessels or berths exceeding their time limits */
    int berth; /**< berths with too many departures on some day */
} Violations;

/** @brief Memory layout of a solution block, shared by all solutions of an instance. */
typedef struct SolutionLayout {
    size_t size; /**< bytes of the block */
    size_t node, pattern, visit_count, visit_node, vessel, voyage, active, berth, day, departure; /**< offsets */
} SolutionLayout;

/** @brief A complete heuristic solution. */
typedef struct Solution {
    const Instance *inst; /**< instance the solution belongs to */
    Cost cost;            /**< penalised objective \f$C(s)\f$ (paper, Section 4.2) */
    bool infeasible;      /**< true if some violation exceeds the tolerance */
    Violations viol;      /**< violation counts of the last feasibility check */

    /* Views into the memory block. Indexing helpers are given below. */
    Node *node;           /**< [num_nodes] */
    int8_t *pattern;      /**< [num_installations] chosen delivery pattern (-1 = not served) */
    int8_t *visit_count;  /**< [num_periods * num_installations] visits of each installation on each day */
    int *visit_node;      /**< [num_periods * num_installations] node visiting the installation on each day */
    Vessel *vessel;       /**< [num_vessels] */
    Voyage *voyage;       /**< [num_vessels * num_periods] */
    int8_t *active;       /**< [num_berths * num_types] active slots of each (berth, type) */
    BerthLog *berth;      /**< [num_berths] */
    BerthDay *day;        /**< [num_berths * num_periods] */
    VoyageKey *departure; /**< [num_berths * num_periods * slots_per_berth] departure queues */

    size_t block_size;    /**< bytes of @c block */
    unsigned char *block; /**< storage of all the arrays above */
} Solution;

/** @brief Computes the memory layout of the solutions of an instance. */
SolutionLayout sol_layout(const Instance *inst);

/**
 * @brief Allocates a solution and initialises it as the empty solution
 *        (no visit served, every voyage empty, every slot inactive).
 *
 * @param inst instance.
 * @return New solution; release it with sol_free().
 */
[[nodiscard]] Solution *sol_create(const Instance *inst);

/** @brief Releases a solution (accepts @c nullptr). */
void sol_free(Solution *s);

/**
 * @brief Copies @p src into @p dst. Both must belong to the same instance.
 *
 * Only the scalar fields and the memory block are copied, so the cost is one
 * @c memcpy of SolutionLayout::size bytes.
 */
void sol_copy(Solution *dst, const Solution *src);

/** @brief Resets @p s to the empty solution (see sol_create()). */
void sol_clear(Solution *s);

/**
 * @name Indexing helpers
 * They accept a const solution because they only compute addresses inside its
 * block; callers that hold a const solution must not write through the result.
 */
///@{
static inline int sol_vessel_of(const Solution *s, VoyageKey k)
{
    return inst_vessel(s->inst, k.berth, k.type, k.slot);
}
static inline Voyage *sol_voyage(const Solution *s, int vessel, int period)
{
    return &s->voyage[vessel * s->inst->num_periods + period];
}
static inline Voyage *sol_voyage_of(const Solution *s, VoyageKey k)
{
    return sol_voyage(s, sol_vessel_of(s, k), k.period);
}
static inline Vessel *sol_vessel(const Solution *s, int berth, int type, int slot)
{
    return &s->vessel[inst_vessel(s->inst, berth, type, slot)];
}
static inline int8_t *sol_active(const Solution *s, int berth, int type)
{
    return &s->active[berth * s->inst->num_types + type];
}
static inline BerthDay *sol_day(const Solution *s, int berth, int period)
{
    return &s->day[berth * s->inst->num_periods + period];
}
static inline VoyageKey *sol_queue(const Solution *s, int berth, int period)
{
    return &s->departure[(berth * s->inst->num_periods + period) * s->inst->slots_per_berth];
}
///@}

/** @brief Time violation of a voyage: time - max_time (negative = slack). */
static inline int voyage_time_slack(const Voyage *r) { return r->time - r->max_time; }

/**
 * @brief Recomputes the penalised objective of @p s.
 *
 * \f$C(s) = \sum_{\text{vessels}} (\text{cost} + \omega^q E^q + \omega^t O^t)
 *  + \sum_{\text{berths}} (\omega^b \cdot \text{excess departures} + \omega^t E^t)\f$.
 *
 * @param s        solution (its cost field is updated).
 * @param w_load   weight of capacity violations \f$\omega^q\f$.
 * @param w_time   weight of time violations \f$\omega^t\f$.
 * @param w_berth  weight of excess departures per berth and day.
 */
void sol_update_cost(Solution *s, float w_load, float w_time, float w_berth);

/**
 * @brief Objective of the paper, Eq. (1): fixed cost of the active vessels plus
 *        the sailing costs, recomputed exactly from the voyages.
 *
 * @return The objective in reported units (internal cost / ::COST_SCALE).
 */
double sol_objective(const Solution *s);

/**
 * @brief Counts the violations of @p s and sets Solution::infeasible.
 *
 * A violation is counted only when it exceeds a tolerance of 10 internal
 * units, which absorbs the rounding of the incremental computations.
 */
void sol_check_feasibility(Solution *s);

/**
 * @brief Removes the gaps left by vessel slots that became inactive.
 *
 * Inactive slots inside the active prefix of a (berth, type) are replaced by
 * the last active slot, so active slots stay contiguous. The time fields of
 * the empty voyages of active vessels are reset to zero.
 */
void sol_compact(Solution *s);

/** @brief Number of active vessels of each type, written to @p count[num_types]. */
void sol_fleet(const Solution *s, int *count);

/**
 * @brief Prints a readable description of the solution: cost, fleet, the
 *        voyages of every vessel and the berth schedules.
 */
void sol_print(const Solution *s, FILE *out);

#endif /* PSVPP_HEURISTIC_SOLUTION_H */
