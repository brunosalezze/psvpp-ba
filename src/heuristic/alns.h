// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file alns.h
 * @brief Adaptive large neighbourhood search with multiple starts and spaced
 *        local searches (paper, Section 4, Algorithms 1 and 2).
 */
#ifndef PSVPP_HEURISTIC_ALNS_H
#define PSVPP_HEURISTIC_ALNS_H

#include "heuristic/operators.h"
#include "heuristic/params.h"
#include "heuristic/search.h"
#include "heuristic/solution.h"

#include <stdio.h>

/** @brief Statistics of one heuristic run. */
typedef struct AlnsReport {
    long iterations;                        /**< ALNS iterations performed */
    long feasible_iterations;               /**< iterations whose candidate was feasible */
    int best_start;                         /**< index of the start that produced the final seed solution */
    double start_cost;                      /**< cost of the best start (before the final ALNS) */
    long removal_uses[NUM_REMOVAL_OPS];     /**< times each removal operator was applied */
    long insertion_uses[NUM_INSERTION_OPS]; /**< times each insertion operator was applied */
    double removal_weight[NUM_REMOVAL_OPS]; /**< final adaptive weights */
    double insertion_weight[NUM_INSERTION_OPS];
} AlnsReport;

/**
 * @brief Runs the ALNS of Algorithm 2.
 *
 * Each iteration destroys and repairs a copy of the current solution with a
 * pair of operators chosen by the adaptive layer (roulette wheel on weights
 * updated every AlnsParams::segment iterations with Eq. 38) and accepts it with
 * the simulated-annealing criterion (Section 4.7). Every
 * AlnsParams::ls_interval iterations the most promising solution (the last
 * local-search result or any cheaper candidate since) is improved by
 * local_search() and becomes the candidate. Penalty
 * weights are adapted every segment (search_update_penalties()).
 *
 * @param S        search context; the penalty weights are used as they are
 *                 (alns_multistart() resets them before every start).
 * @param io       input: initial solution, with its cost computed under the
 *                 current weights; output: best feasible solution
 *                 found, or the best penalised solution if none was feasible.
 * @param max_iter number of iterations \f$N\f$.
 * @param report   statistics are accumulated here (may be @c nullptr).
 */
void alns_run(Search *S, Solution *io, long max_iter, AlnsReport *report);

/**
 * @brief Multi-start ALNS of Algorithm 1.
 *
 * Builds AlnsParams::starts random solutions, improves each with
 * AlnsParams::start_iterations ALNS iterations, then runs the ALNS on the best
 * one for the remaining AlnsParams::iterations - starts * start_iterations
 * iterations. Feasible solutions are always preferred to infeasible ones.
 *
 * @param S      search context.
 * @param best   receives the best solution found.
 * @param report receives the run statistics (may be @c nullptr).
 * @param log    progress messages are written here when not @c nullptr.
 */
void alns_multistart(Search *S, Solution *best, AlnsReport *report, FILE *log);

#endif /* PSVPP_HEURISTIC_ALNS_H */
