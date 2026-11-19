// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file exact.h
 * @brief Public entry points of the exact methods (paper, Sections 2 and 3).
 *
 * The voyage enumeration (`routes`) needs no external library. The
 * branch-and-cut (`bc`) and the complete model (`complete`) need Gurobi; when
 * the program is built without it (`make GUROBI=0`) these commands report an
 * error.
 */
#ifndef PSVPP_EXACT_EXACT_H
#define PSVPP_EXACT_EXACT_H

#include "common/instance.h"

/** Largest voyage size accepted by --max-route-size (the TSP is solved by dynamic programming). */
enum { MAX_ROUTE_SIZE_LIMIT = 15 };

/** @brief Parameters of the exact methods. */
typedef struct ExactParams {
    double time_limit;      /**< overall time limit of the MIP solver (s); the paper used 86,400 */
    double gap;             /**< relative gap at which the solver stops (1e-4, Gurobi's default; the original
                                 code stopped the branch-and-cut at 1e-3) */
    int threads;            /**< Gurobi threads of the main model (reduced or complete), 0 = all cores */
    int subproblem_threads; /**< bc: Gurobi threads of the models that answer when a scheduling algorithm
                                 gives up (see separation.h), 0 = all cores; never more than @c threads */
    int max_route_size;     /**< maximum installations per voyage (8, as in the paper) */
    double min_usage;       /**< bc only: on non-clustered instances drop voyages using less than this
                                 percentage of the vessel capacity (0 = exact; see docs/algorithms.md) */
    bool solver_log;        /**< print the Gurobi log */
    bool root_only;         /**< bc: stop after the LP relaxation of the reduced model (a lower bound only) */
} ExactParams;

/** @brief Default parameters. */
ExactParams exact_default_params(void);

/** @brief `routes` command: enumerates the voyages and prints statistics. */
int exact_routes_command(const Instance *inst, const ExactParams *par, bool quiet);

/** @brief `bc` command: branch-and-cut on the reduced model. */
int exact_bc_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output);

/** @brief `complete` command: complete model with overlap cuts. */
int exact_complete_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output);

#endif /* PSVPP_EXACT_EXACT_H */
