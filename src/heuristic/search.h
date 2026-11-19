// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file search.h
 * @brief State shared by all heuristic components during one run.
 *
 * A ::Search bundles everything that the original implementation kept in
 * global variables: the random number generator, the adaptive penalty weights
 * and their feasibility statistics, and scratch buffers sized for the
 * instance. One ::Search per run makes runs independent and thread-safe.
 */
#ifndef PSVPP_HEURISTIC_SEARCH_H
#define PSVPP_HEURISTIC_SEARCH_H

#include "common/instance.h"
#include "common/rng.h"
#include "heuristic/params.h"
#include "heuristic/solution.h"

#include <stdint.h>

/** @brief Change of the waiting time of a queued voyage caused by a pending move. */
typedef struct QueueShift {
    int8_t period; /**< departure day of the affected voyage */
    int8_t type;   /**< vessel type of the affected voyage */
    int8_t slot;   /**< vessel slot of the affected voyage */
    int d_wait;    /**< change of its waiting time (s) */
    int d_excess;  /**< change of its vessel's time violation (s) */
} QueueShift;

/** @brief Candidate insertion position and its cost (see insertion.h). */
typedef struct InsertOption {
    int before; /**< node after which the visit would be inserted */
    Cost cost;  /**< penalised insertion cost (or regret value) */
} InsertOption;

/** @brief Penalty weights \f$\omega\f$ of the three kinds of violations. */
typedef struct PenaltyWeights {
    float load;  /**< \f$\omega^q\f$: per load unit beyond the vessel capacity (and see Search::size_excess) */
    float time;  /**< \f$\omega^t\f$: per second of time violation */
    float berth; /**< per departure beyond ::MAX_DEPARTURES_PER_BERTH on a berth-day */
} PenaltyWeights;

/** @brief Feasibility statistics of the current segment of iterations. */
typedef struct FeasibilityStats {
    int total; /**< solutions observed */
    int load;  /**< ... with capacity violations */
    int time;  /**< ... with time violations */
    int berth; /**< ... with too many departures */
} FeasibilityStats;

/** @brief Run-wide state of the heuristic. */
typedef struct Search {
    const Instance *inst;   /**< instance being solved */
    const AlnsParams *par;  /**< parameters */
    Rng rng;                /**< random stream of this run */
    PenaltyWeights w;       /**< current penalty weights */
    int size_excess;        /**< load counted as a capacity violation for every installation of a voyage
                                 beyond AlnsParams::max_route_size (weighed by PenaltyWeights::load) */
    FeasibilityStats stats; /**< statistics since the last weight update */

    /* scratch buffers of the move evaluator (moves.c) */
    int *shift;        /**< [num_periods] change of the carried berth time per day */
    int *d_sink;       /**< [num_periods] change of the absorbed time per day */
    bool shift_dirty;  /**< shift or d_sink hold non-zero entries (of the last evaluated move) */
    QueueShift *queue; /**< [num_periods * slots_per_berth] waiting-time changes */

    /* scratch buffers of the insertion evaluator (insertion.c) */
    int *cand;              /**< [num_installations + num_vessels] candidate positions */
    InsertOption *day_best; /**< [num_periods] best option of each day */
    Cost *pattern_cost;     /**< [max number of patterns] cost of each delivery pattern */
    int *choice_buf;        /**< [2 * num_periods] storage for two ::PatternChoice position arrays */

    /* scratch buffers of the construction (construction.c) */
    int *tour;      /**< [num_periods * num_installations] giant tour of each day */
    int *tour_size; /**< [num_periods] visits in each giant tour */

    /* scratch buffers of the local search (local_search.c) */
    Solution *undo; /**< snapshot restored when a tentative move is rejected */
    int *ls_cand;   /**< [num_installations + num_vessels] neighbour positions */
    int *ls_order;  /**< [max(num_visits, num_installations)] random visiting order */
} Search;

/**
 * @brief Initialises a search context.
 *
 * @param S     context to initialise.
 * @param inst  instance (must outlive the context).
 * @param par   parameters (must outlive the context).
 * @param seed  seed of the random stream.
 */
void search_init(Search *S, const Instance *inst, const AlnsParams *par, uint64_t seed);

/** @brief Releases the buffers of a search context. */
void search_free(Search *S);

/** @brief Restores the initial penalty weights and clears the statistics. */
void search_reset_penalties(Search *S);

/** @brief Records the violations of @p s in the feasibility statistics. */
void search_observe(Search *S, const Solution *s);

/**
 * @brief Adapts the penalty weights to the statistics of the last segment.
 *
 * For each kind of violation, the weight is multiplied by
 * AlnsParams::penalty_increase when more than target + 5% of the observed
 * solutions violated it, and by AlnsParams::penalty_decrease when fewer than
 * target - 5% did, within the configured bounds. The statistics are cleared.
 */
void search_update_penalties(Search *S);

/** @brief Recomputes the penalised cost of @p s with the current weights. */
static inline void search_update_cost(const Search *S, Solution *s)
{
    sol_update_cost(s, S->w.load, S->w.time, S->w.berth);
}

#endif /* PSVPP_HEURISTIC_SEARCH_H */
