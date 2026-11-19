// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file separation.h
 * @brief Feasibility checks of an integer solution of the reduced model
 *        (paper, Section 3): the building blocks of the separation procedure
 *        (Algorithm 1) and of Reduce (Algorithm 2).
 *
 * An incumbent of the reduced model fixes, for every selected voyage, its
 * departure day and vessel type, plus the fleet size of each type. It is
 * feasible for the complete model if every voyage can also receive a vessel,
 * a berth and a departing position such that no vessel is late for its next
 * departure. The checks below decide this, cheapest first:
 *  1. sep_berth_relaxation(): the berths alone (vessel timing ignored);
 *  2. sep_colouring(): the vessels alone (berths ignored) — a minimum vertex
 *     colouring of the circular-arc graph of the voyages;
 *  3. sep_subproblem(): both together (*SubProblem* of Section 3.1).
 *
 * The paper solves these problems as MIP models; here they are solved by the
 * combinatorial algorithms of schedule.h, which need no solver. An algorithm
 * that reaches its work limit without an answer hands the plan over to
 * ExactContext::fallback when one is installed.
 */
#ifndef PSVPP_EXACT_SEPARATION_H
#define PSVPP_EXACT_SEPARATION_H

#include "exact/exact.h"
#include "exact/plan.h"
#include "exact/routes.h"
#include "exact/schedule.h"

/**
 * @brief Exact answers for the plans on which an algorithm reaches its work limit (optional).
 *
 * The branch-and-cut installs the MIP models of subproblem_models.h here.
 */
typedef struct ExactFallback {
    bool (*berth_relaxation)(void *context, Plan *plan);                /**< as sep_berth_relaxation() */
    int (*subproblem)(void *context, Plan *plan, bool keep_assignment); /**< 1 feasible, 0 infeasible, -1 unknown */
    void *context;
} ExactFallback;

/** @brief Shared state of the checks. Create with exact_context_init(), release with exact_context_free(). */
typedef struct ExactContext {
    const Instance *inst;
    const ExactParams *par;
    const RoutePool *pool;
    Scheduler *scheduler; /**< working memory of the scheduling algorithms */
    ExactFallback fallback;
    int *work; /**< scratch arrays of the checks */
    bool *include;
    int work_capacity;
    /* statistics */
    long incumbents;           /**< integer solutions checked */
    long berth_cuts;           /**< cuts from the berth relaxation */
    long colouring_cuts;       /**< cuts (33) from the fleet colouring */
    long subproblem_cuts;      /**< cuts (32) from the feasibility subproblem */
    long subproblems;          /**< subproblems solved (including those of Reduce) */
    long fallbacks;            /**< checks answered by the fallback after an algorithm gave up */
    long unproven_subproblems; /**< subproblems left without answer (counted as infeasible) */
    double check_seconds;      /**< wall time spent in the checks */
} ExactContext;

/** @brief Prepares a context for the checks of plans over @p pool. */
void exact_context_init(ExactContext *X, const Instance *inst, const ExactParams *par, const RoutePool *pool);

/** @brief Releases the memory of a context. */
void exact_context_free(ExactContext *X);

/**
 * @brief Berth-only relaxation of the subproblem.
 *
 * Looks for a berth and a departing position for every voyage such that no
 * berth carries more than TOL - 1 into any day, whatever the waiting times of
 * the vessels. When there is none, Trip::conflict marks a minimal set of
 * voyages that cannot be scheduled on their own: whole days, then single
 * voyages, are removed from the plan as long as it stays infeasible.
 *
 * @return true if the berth schedule is feasible.
 */
bool sep_berth_relaxation(ExactContext *X, Plan *plan);

/**
 * @brief Vessel-only relaxation for one vessel type: can the voyages of type
 *        @p type be performed by Plan::fleet[@p type] vessels without overlaps?
 *
 * When they cannot, @p needed receives the smallest fleet that can (the
 * chromatic number of their circular-arc graph) and Trip::conflict marks,
 * among the voyages of the type, a minimal set that already needs that fleet.
 *
 * @return true if the fleet is sufficient.
 */
bool sep_colouring(ExactContext *X, Plan *plan, int type, int *needed);

/**
 * @brief Feasibility subproblem (paper, Section 3.1), solved by schedule_plan().
 *
 * @param X               context.
 * @param plan            voyages and fleet; on success with @p keep_assignment
 *                        the vessel, berth, position and waiting times are stored.
 * @param keep_assignment store the assignment found.
 * @return true if the plan is feasible for the complete model.
 */
bool sep_subproblem(ExactContext *X, Plan *plan, bool keep_assignment);

/**
 * @brief Reduce (Algorithm 2): shrinks an infeasible plan to a subset of
 *        voyages that is still infeasible for the subproblem, minimal unless
 *        some of its subproblems could not be answered (their voyages are kept).
 */
void sep_reduce(ExactContext *X, Plan *plan);

/**
 * @brief Verifies the stored assignment of a plan directly from the rules of
 *        the complete model (schedule_check()).
 */
bool sep_verify(ExactContext *X, const Plan *plan);

#endif /* PSVPP_EXACT_SEPARATION_H */
