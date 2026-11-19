// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file insertion.h
 * @brief Best insertion of an installation over its delivery patterns, the
 *        operator \f$\Phi(u, s)\f$ of the paper (Section 4.4).
 *
 * Inserting an installation means choosing one of its delivery patterns and,
 * for every day of the pattern, a position in some voyage of that day. The
 * candidate positions of a day are the visits of the installation's nearest
 * neighbours served that day, plus the head node of every active vessel and of
 * the first inactive slot of each (berth, type), so a new voyage or a new
 * vessel can always be opened. Days are evaluated independently.
 */
#ifndef PSVPP_HEURISTIC_INSERTION_H
#define PSVPP_HEURISTIC_INSERTION_H

#include "heuristic/search.h"
#include "heuristic/solution.h"

/** @brief An installation together with the chosen pattern and positions. */
typedef struct PatternChoice {
    int inst;    /**< installation */
    int pattern; /**< chosen delivery pattern */
    Cost value;  /**< total insertion cost of the pattern (or its regret value) */
    int *before; /**< [num_periods] node after which the visit of each day is inserted */
} PatternChoice;

/**
 * @brief Collects the candidate insertion positions of an installation on one day.
 *
 * @param s      solution.
 * @param inst   installation to insert.
 * @param period day.
 * @param out    receives the candidate nodes (capacity: num_installations + num_vessels).
 * @return Number of candidates written to @p out.
 */
int insertion_candidates(const Solution *s, int inst, int period, int *out);

/**
 * @brief Finds the best delivery pattern and positions for an unrouted installation.
 *
 * With @p regret_k < 2 the value of a day is the cheapest insertion cost and
 * the pattern with the cheapest total is chosen (greedy insertion). With
 * @p regret_k >= 2 the value of a day is its regret,
 * \f$\sum_{i=2}^{k} (c_1 - c_i)\f$ over its @p regret_k cheapest positions
 * (a non-positive number), and the pattern with the largest total regret is
 * chosen (k-regret insertion, IO3). In both cases the positions are the
 * cheapest ones of each day.
 *
 * @param S        search context.
 * @param s        solution (not modified).
 * @param choice   input: PatternChoice::inst; output: the other fields.
 * @param regret_k regret level (see above).
 */
void insertion_best_pattern(Search *S, Solution *s, PatternChoice *choice, int regret_k);

/**
 * @brief Inserts an installation according to a ::PatternChoice.
 *
 * @return Change of the penalised objective.
 */
Cost insertion_apply(Search *S, Solution *s, const PatternChoice *choice);

/**
 * @brief Removes every visit of an installation from the solution.
 *
 * @return Change of the penalised objective.
 */
Cost insertion_remove_installation(Search *S, Solution *s, int inst);

/**
 * @brief Cost change of removing every visit of an installation, each visit
 *        evaluated independently (the solution is not modified).
 */
Cost insertion_removal_cost(Search *S, Solution *s, int inst);

#endif /* PSVPP_HEURISTIC_INSERTION_H */
