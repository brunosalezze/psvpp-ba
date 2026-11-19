// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file exact.c
 * @brief Commands of the exact methods (see exact.h).
 */
#include "exact/exact.h"

#include "exact/routes.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef PSVPP_WITH_GUROBI
#include "exact/models.h"
#endif

ExactParams exact_default_params(void)
{
    return (ExactParams){
        .time_limit = 86400.0,
        .gap = 1e-4,
        .threads = 0,
        .subproblem_threads = 1,
        .max_route_size = 8,
        .min_usage = 0.0,
        .solver_log = false,
        .root_only = false,
    };
}

int exact_routes_command(const Instance *inst, const ExactParams *par, bool quiet)
{
    if (!quiet)
        inst_print_summary(inst);
    RoutePool *pool = route_pool_build(inst, par->max_route_size);
    route_pool_print_stats(pool);
    route_pool_free(pool);
    return EXIT_SUCCESS;
}

#ifdef PSVPP_WITH_GUROBI

int exact_bc_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output)
{
    return solve_reduced_model(inst, par, quiet, output);
}

int exact_complete_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output)
{
    return solve_complete_model(inst, par, quiet, output);
}

#else

/* Reports that the exact models are unavailable in this build; returns EXIT_FAILURE. */
static int no_gurobi(void)
{
    fprintf(stderr, "error: this build has no Gurobi support; rebuild with Gurobi installed "
                    "(see README.md, section 'Building')\n");
    return EXIT_FAILURE;
}

int exact_bc_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output)
{
    (void)inst, (void)par, (void)quiet, (void)output;
    return no_gurobi();
}

int exact_complete_command(const Instance *inst, const ExactParams *par, bool quiet, const char *output)
{
    (void)inst, (void)par, (void)quiet, (void)output;
    return no_gurobi();
}

#endif
