// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file local_search.h
 * @brief Local search of the ALNS (paper, Section 4.5 and Algorithm 3).
 *
 * Two first-improvement neighbourhoods are alternated, starting with a random one:
 *  - **VI (voyage improvement)**: for a visit @c u and a candidate position
 *    @c v of the same day, exchange the sequence of up to 3 visits starting at
 *    @c u with the sequence of up to 2 visits following @c v (including the
 *    pure relocation of @c u), also trying the reversed sequences;
 *  - **PI (pattern improvement)**: remove an installation and reinsert it with
 *    the best pattern and positions (insertion_best_pattern()).
 * A move is accepted when it reduces the penalised cost by more than a small
 * tolerance (::LS_IMPROVEMENT_TOLERANCE internal cost units).
 */
#ifndef PSVPP_HEURISTIC_LOCAL_SEARCH_H
#define PSVPP_HEURISTIC_LOCAL_SEARCH_H

#include "heuristic/search.h"
#include "heuristic/solution.h"

/** Minimum cost decrease for a local-search move to count as an improvement. */
enum { LS_IMPROVEMENT_TOLERANCE = 100 };

/**
 * @brief One pass of the voyage-improvement neighbourhood (VI).
 *
 * Visits are examined in random order; for each, the candidate positions of
 * its day are examined in random order. The first improving exchange is kept.
 *
 * @return true if an improving move was applied.
 */
bool ls_voyage_improvement(Search *S, Solution *s);

/**
 * @brief One pass of the pattern-improvement neighbourhood (PI).
 *
 * Installations are examined in random order; the first one whose best
 * reinsertion improves the solution is kept reinserted, otherwise the original
 * positions are restored.
 *
 * @return true if an improving move was applied.
 */
bool ls_pattern_improvement(Search *S, Solution *s);

/**
 * @brief Local search of Algorithm 3: alternates VI and PI until neither
 *        improves (or AlnsParams::ls_max_rounds steps are done), then
 *        recomputes the penalised cost of @p s.
 */
void local_search(Search *S, Solution *s);

#endif /* PSVPP_HEURISTIC_LOCAL_SEARCH_H */
