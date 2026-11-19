// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file models.h
 * @brief The two MIP formulations of the paper (built with Gurobi).
 */
#ifndef PSVPP_EXACT_MODELS_H
#define PSVPP_EXACT_MODELS_H

#include "common/instance.h"
#include "exact/exact.h"

/**
 * @brief Branch-and-cut on the reduced model (paper, Sections 2.2 and 3).
 *
 * Variables: \f$x_{rtv}\f$ (voyage r departs on day t with a vessel of type
 * v), \f$y_v\f$ (fleet of type v), \f$z_{ip}\f$ (pattern p of installation
 * i) and \f$c_{vtl}\f$ (some voyage of type v departs on day t and spans l
 * days). Constraints (24)-(27) are stated explicitly; every integer solution
 * is checked in a lazy-constraint callback by the separation procedure of
 * Algorithm 1, which adds cuts (32) and (33)-(34) when it is not feasible for the
 * complete model.
 *
 * The solver is run in rounds of increasing time limits; between rounds the
 * cuts found so far are added to the model as ordinary constraints.
 *
 * @return @c EXIT_SUCCESS, or @c EXIT_FAILURE if no solution was found.
 */
int solve_reduced_model(const Instance *inst, const ExactParams *par, bool quiet, const char *output);

/**
 * @brief Complete model (paper, Section 2.1) with the overlap cuts (22)
 *        added in a lazy-constraint callback.
 *
 * @return @c EXIT_SUCCESS, or @c EXIT_FAILURE if no solution was found.
 */
int solve_complete_model(const Instance *inst, const ExactParams *par, bool quiet, const char *output);

#endif /* PSVPP_EXACT_MODELS_H */
