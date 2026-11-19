// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file subproblem_models.h
 * @brief The berth relaxation and the feasibility subproblem as MIP models
 *        (paper, Section 3.1), solved with Gurobi.
 *
 * These are the models of the paper. The separation normally answers the same
 * questions with the algorithms of schedule.h, which are about a hundred times
 * faster; the models are called when an algorithm reaches its work limit
 * without an answer (dense plans whose infeasibility is hard to prove by
 * search), so that the branch-and-cut stays exact.
 */
#ifndef PSVPP_EXACT_SUBPROBLEM_MODELS_H
#define PSVPP_EXACT_SUBPROBLEM_MODELS_H

#include "common/rng.h"
#include "exact/grb.h"
#include "exact/plan.h"
#include "exact/routes.h"

/** @brief State shared by the models. */
typedef struct ModelContext {
    const Instance *inst;
    const RoutePool *pool;
    GRBenv *env; /**< Gurobi environment */
    int threads; /**< Gurobi threads of each model */
    Rng rng;     /**< random order in which a berth conflict is reduced */
} ModelContext;

/**
 * @brief Berth-only relaxation: assigns the voyages to berths and positions
 *        minimising the largest waiting time carried into a day beyond TOL - 1.
 *
 * When the minimum is positive, Trip::conflict marks a set of voyages that
 * cannot be scheduled on their own: the voyages of the most occupied day and
 * of the preceding days that carried waiting time into it when these are
 * infeasible by themselves (checked with a second model), every voyage of the
 * plan otherwise.
 *
 * @return true if the berth schedule is feasible.
 */
bool model_berth_relaxation(ModelContext *X, Plan *plan);

/**
 * @brief Feasibility subproblem, solved by decomposition: a vessel assignment
 *        (AV) is generated, the slack of every voyage until its vessel's next
 *        departure is derived, and a berth assignment (AB) respecting these
 *        slacks is sought. When AB fails, a minimal set of conflicting voyages
 *        is identified and the vessel assignment of these voyages and of the
 *        voyages that follow them on their vessels (which set the slacks) is
 *        forbidden in AV (no-good cut); conflict sets are cached in a tabu list.
 *
 * @param X               context.
 * @param plan            voyages and fleet; on success with @p keep_assignment
 *                        the vessel, berth, position and waiting times are stored.
 * @param keep_assignment store the assignment found.
 * @return 1 if the plan is feasible for the complete model, 0 if it is not,
 *         -1 if the iteration limit was reached without an answer.
 */
int model_subproblem(ModelContext *X, Plan *plan, bool keep_assignment);

#endif /* PSVPP_EXACT_SUBPROBLEM_MODELS_H */
