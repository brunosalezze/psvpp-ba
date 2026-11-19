// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file alns.c
 * @brief ALNS main loop, adaptive layer and multi-start (see alns.h).
 */
#include "heuristic/alns.h"

#include "common/util.h"
#include "heuristic/construction.h"
#include "heuristic/local_search.h"

#include <math.h>
#include <string.h>

AlnsParams alns_default_params(void)
{
    return (AlnsParams){
        .starts = 50,
        .start_iterations = 100,
        .iterations = 300000,
        .ls_interval = 50,
        .segment = 50,
        .t0 = 1.0 / 50.0,
        .tf = 1.0 / 500.0,
        .removal_min = 0.10,
        .removal_max = 0.25,
        .removal_extra = 3,
        .determinism = 4.0,
        .insertion_sample = 5,
        .adaptive = true,
        .reaction = 0.05,
        .score_best = 5.0,
        .score_better = 0.0,
        .score_accept = 10.0,
        .penalty_start = 700.0,
        .penalty_increase = 1.15,
        .penalty_decrease = 0.90,
        .target_infeasible = 0.60,
        .penalty_min = 10.0,
        .penalty_max_load = 2500.0,
        .penalty_max_time = 1500.0,
        .penalty_max_berth = 2000000.0,
        .max_route_size = 8,
        .ls_max_rounds = 80,
    };
}

/* ------------------------------------------------------------------------- */
/*  Adaptive layer (paper, Section 4.6)                                      */
/* ------------------------------------------------------------------------- */

enum { MAX_OPERATORS = 8 };

/** @brief Roulette wheel over a family of operators. */
typedef struct Roulette {
    int n;                        /**< number of operators */
    double weight[MAX_OPERATORS]; /**< selection weights \f$\phi_i\f$ */
    double score[MAX_OPERATORS];  /**< scores \f$\pi_i\f$ collected in the current segment */
    int used[MAX_OPERATORS];      /**< uses \f$\xi_i\f$ in the current segment */
} Roulette;

/* Initialises a roulette wheel over @p n operators, all with weight 1. */
static void roulette_init(Roulette *r, int n)
{
    *r = (Roulette){.n = n};
    for (int i = 0; i < n; i++)
        r->weight[i] = 1.0;
}

/* Draws an operator: proportionally to the weights when the adaptive layer is on, uniformly otherwise.
 * Returns the operator index. */
static int roulette_pick(Roulette *r, Search *S)
{
    if (!S->par->adaptive)
        return rng_int(&S->rng, 0, r->n - 1);
    double total = 0;
    for (int i = 0; i < r->n; i++)
        total += r->weight[i];
    if (!(total > 0)) /* every weight decayed to zero (possible with zero scores): choose uniformly */
        return rng_int(&S->rng, 0, r->n - 1);
    double x = rng_unit(&S->rng) * total;
    for (int i = 0; i < r->n - 1; i++) {
        x -= r->weight[i];
        if (x < 0)
            return i;
    }
    return r->n - 1;
}

/* Eq. (38): phi_i = (1 - tau) phi_i + tau * pi_i / xi_i for the operators used in the segment. */
static void roulette_update(Roulette *r, double reaction)
{
    for (int i = 0; i < r->n; i++) {
        if (r->used[i] > 0)
            r->weight[i] = (1 - reaction) * r->weight[i] + reaction * r->score[i] / r->used[i];
        r->score[i] = 0;
        r->used[i] = 0;
    }
}

/* ------------------------------------------------------------------------- */
/*  ALNS (Algorithm 2)                                                       */
/* ------------------------------------------------------------------------- */

/** Outcome of the acceptance test of one iteration. */
typedef enum Outcome { NEW_BEST, IMPROVED, ACCEPTED, REJECTED } Outcome;

void alns_run(Search *S, Solution *io, long max_iter, AlnsReport *report)
{
    const Instance *in = S->inst;
    const AlnsParams *par = S->par;
    Solution *cur = sol_create(in), *cand = sol_create(in), *best = sol_create(in), *prom = sol_create(in);
    InstallationPool pool;
    pool_init(&pool, in);
    Roulette rem, ins;
    roulette_init(&rem, NUM_REMOVAL_OPS);
    roulette_init(&ins, NUM_INSERTION_OPS);

    /* Temperature schedule (Eqs. 39-41): theta_0 = C(s) theta'_0, cooled to C(s) theta'_F. */
    const float t0 = (float)par->t0, tf = (float)par->tf;
    float temp = io->cost * t0;
    const float cooling = (float)pow(tf / t0, 1 / ((float)max_iter - 1.0f));

    sol_copy(best, io);
    sol_copy(cur, best);
    /* prom is the "most promising" solution s_prom: the last local-search result, replaced by any
     * cheaper candidate. The paper's Algorithm 2 resets C(s_prom) after each local search; the
     * original implementation effectively kept it (its cost was recomputed right after the reset),
     * and keeping it gives clearly better solutions, see docs/changes-from-original.md. */
    bool have_prom = false;
    for (long it = 1; it <= max_iter; it++) {
        sol_copy(cand, cur);
        const int r = roulette_pick(&rem, S);
        apply_removal(S, cand, &pool, (RemovalOp)r);
        const int k = roulette_pick(&ins, S);
        apply_insertion(S, cand, &pool, (InsertionOp)k);
        search_update_cost(S, cand);

        /* Spaced local search on the most promising solution. */
        if (!have_prom || prom->cost > cand->cost) {
            sol_copy(prom, cand);
            have_prom = true;
        }
        if (it % par->ls_interval == 0) {
            local_search(S, prom);
            sol_copy(cand, prom);
        }
        sol_compact(cand);
        sol_check_feasibility(cand);

        /* Penalty weights follow the feasibility of the candidates and of the best solution. */
        search_observe(S, cand);
        search_observe(S, best);
        if (it % par->segment == 0) {
            search_update_penalties(S);
            if (have_prom)
                search_update_cost(S, prom);
            search_update_cost(S, cand);
            search_update_cost(S, cur);
            search_update_cost(S, best);
        }

        if (!cand->infeasible) {
            if (report)
                report->feasible_iterations++;
            if (cand->cost <= io->cost || io->infeasible)
                sol_copy(io, cand);
        }

        Outcome outcome;
        if (best->cost >= cand->cost) {
            sol_copy(cur, cand);
            sol_copy(best, cand);
            outcome = NEW_BEST;
        } else if (cur->cost >= cand->cost) {
            sol_copy(cur, cand);
            outcome = IMPROVED;
        } else if (rng_unit(&S->rng) < expf((cur->cost - cand->cost) / temp)) {
            sol_copy(cur, cand);
            outcome = ACCEPTED;
        } else {
            outcome = REJECTED;
        }

        const double score = outcome == NEW_BEST   ? par->score_best
                             : outcome == IMPROVED ? par->score_better
                             : outcome == ACCEPTED ? par->score_accept
                                                   : 0.0;
        rem.score[r] += score;
        ins.score[k] += score;
        rem.used[r]++;
        ins.used[k]++;
        if (it % par->segment == 0) {
            roulette_update(&rem, par->reaction);
            roulette_update(&ins, par->reaction);
        }
        if (report) {
            report->iterations++;
            report->removal_uses[r]++;
            report->insertion_uses[k]++;
        }
        temp = cooling * temp;
    }
    if (io->infeasible)
        sol_copy(io, best);

    if (report) {
        for (int i = 0; i < NUM_REMOVAL_OPS; i++)
            report->removal_weight[i] = rem.weight[i];
        for (int i = 0; i < NUM_INSERTION_OPS; i++)
            report->insertion_weight[i] = ins.weight[i];
    }
    pool_free(&pool);
    sol_free(cur);
    sol_free(cand);
    sol_free(best);
    sol_free(prom);
}

/* ------------------------------------------------------------------------- */
/*  Multi-start (Algorithm 1)                                                */
/* ------------------------------------------------------------------------- */

/* Feasible solutions first, then by penalised cost. */
static bool better(const Solution *a, const Solution *b)
{
    if (a->infeasible != b->infeasible)
        return !a->infeasible;
    return a->cost <= b->cost;
}

void alns_multistart(Search *S, Solution *best, AlnsReport *report, FILE *log)
{
    const AlnsParams *par = S->par;
    const long final_iter = (long)par->iterations - (long)par->starts * par->start_iterations;
    if (final_iter < 0)
        die("the total number of iterations (%d) is smaller than starts x iterations per start (%d x %d)",
            par->iterations, par->starts, par->start_iterations);
    if (report)
        *report = (AlnsReport){.best_start = -1};

    if (par->starts == 0) {
        search_reset_penalties(S);
        construct_random_solution(S, best);
    } else {
        Solution *trial = sol_create(S->inst);
        for (int k = 0; k < par->starts; k++) {
            search_reset_penalties(S);
            construct_random_solution(S, trial);
            alns_run(S, trial, par->start_iterations, report);
            if (k == 0 || better(trial, best)) {
                sol_copy(best, trial);
                if (report)
                    report->best_start = k;
            }
            if (log)
                fprintf(log, "start %3d/%d: %s%.6f (best %.6f)\n", k + 1, par->starts, trial->infeasible ? "*" : "",
                        trial->cost / COST_SCALE, best->cost / COST_SCALE);
        }
        sol_free(trial);
    }
    if (report)
        report->start_cost = best->cost;

    if (final_iter > 0) {
        search_reset_penalties(S);
        search_update_cost(S, best);
        alns_run(S, best, final_iter, report);
    }
    if (log)
        fprintf(log, "final: %s%.6f\n", best->infeasible ? "*" : "", best->cost / COST_SCALE);
}
