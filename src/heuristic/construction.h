// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file construction.h
 * @brief Random initial solutions: the Fix and Split phases (paper, Section 4.3).
 */
#ifndef PSVPP_HEURISTIC_CONSTRUCTION_H
#define PSVPP_HEURISTIC_CONSTRUCTION_H

#include "heuristic/search.h"
#include "heuristic/solution.h"

/**
 * @brief Builds a random initial solution.
 *
 * - **Fix**: every installation receives a random delivery pattern; the visits
 *   of each day form a randomly shuffled giant tour.
 * - **Split**: each giant tour is cut, in order, into voyages of randomly
 *   chosen (berth, vessel type) slots, filling each vessel up to a fraction of
 *   its capacity. The fraction starts at 75% and grows by 25% whenever the
 *   available slots are not enough to serve every visit.
 *
 * The result may violate time constraints (they are penalised), never the
 * structure. Its feasibility flags and penalised cost are up to date.
 *
 * @param S search context (random stream and penalty weights).
 * @param s solution to overwrite.
 */
void construct_random_solution(Search *S, Solution *s);

#endif /* PSVPP_HEURISTIC_CONSTRUCTION_H */
