// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file search.c
 * @brief Search context and adaptive feasibility penalties (see search.h).
 */
#include "heuristic/search.h"

#include "common/util.h"

#include <stdlib.h>

void search_init(Search *S, const Instance *inst, const AlnsParams *par, uint64_t seed)
{
    const int cand = inst->num_installations + inst->num_vessels + 1;
    *S = (Search){.inst = inst, .par = par};
    rng_seed(&S->rng, seed, 0);
    S->shift = XCALLOC_ARRAY(int, inst->num_periods);
    S->d_sink = XCALLOC_ARRAY(int, inst->num_periods);
    S->queue = XCALLOC_ARRAY(QueueShift, inst->num_periods * inst->slots_per_berth);
    S->cand = XCALLOC_ARRAY(int, cand);
    int max_patterns = 1;
    for (int i = 0; i < inst->num_installations; i++)
        max_patterns = imax(max_patterns, inst->inst[i].num_patterns);
    S->day_best = XCALLOC_ARRAY(InsertOption, inst->num_periods);
    S->pattern_cost = XCALLOC_ARRAY(Cost, max_patterns);
    S->choice_buf = XCALLOC_ARRAY(int, 2 * inst->num_periods);
    S->tour = XCALLOC_ARRAY(int, inst->num_periods * inst->num_installations);
    S->tour_size = XCALLOC_ARRAY(int, inst->num_periods);
    S->undo = sol_create(inst);
    S->ls_cand = XCALLOC_ARRAY(int, cand);
    S->ls_order = XCALLOC_ARRAY(int, imax(inst->num_visits, inst->num_installations));
    /* An installation beyond the largest voyage size weighs like the average visit of a full voyage
     * of that size on the largest vessel. */
    int max_capacity = 0;
    for (int k = 0; k < inst->num_types; k++)
        max_capacity = imax(max_capacity, inst->type[k].capacity);
    S->size_excess = imax(max_capacity / imax(par->max_route_size, 1), 1000); /* far above the feasibility tolerance */
    search_reset_penalties(S);
}

void search_free(Search *S)
{
    free(S->shift);
    free(S->d_sink);
    free(S->queue);
    free(S->cand);
    free(S->day_best);
    free(S->pattern_cost);
    free(S->choice_buf);
    free(S->tour);
    free(S->tour_size);
    sol_free(S->undo);
    free(S->ls_cand);
    free(S->ls_order);
    *S = (Search){0};
}

void search_reset_penalties(Search *S)
{
    const float w0 = (float)S->par->penalty_start;
    S->w = (PenaltyWeights){.load = w0, .time = w0, .berth = w0};
    S->stats = (FeasibilityStats){0};
}

void search_observe(Search *S, const Solution *s)
{
    S->stats.total++;
    S->stats.load += s->viol.load > 0;
    S->stats.time += s->viol.time > 0;
    S->stats.berth += s->viol.berth > 0;
}

/* Multiplies one weight up or down according to the observed infeasibility rate. */
static void adapt_weight(const AlnsParams *par, int total, int infeasible, float *w, double max_weight)
{
    const float up = (float)par->penalty_increase, down = (float)par->penalty_decrease;
    if (total * (par->target_infeasible + 0.05) < infeasible && *w < max_weight)
        *w *= up;
    if (total * (par->target_infeasible - 0.05) > infeasible && *w > par->penalty_min)
        *w *= down;
}

void search_update_penalties(Search *S)
{
    const AlnsParams *par = S->par;
    adapt_weight(par, S->stats.total, S->stats.load, &S->w.load, par->penalty_max_load);
    adapt_weight(par, S->stats.total, S->stats.time, &S->w.time, par->penalty_max_time);
    adapt_weight(par, S->stats.total, S->stats.berth, &S->w.berth, par->penalty_max_berth);
    S->stats = (FeasibilityStats){0};
}
