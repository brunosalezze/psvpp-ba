// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file params.h
 * @brief Parameters of the multi-start ALNS heuristic (paper, Table 3).
 *
 * The defaults returned by alns_default_params() are the best values reported
 * by the paper's tuning (Section 5.2), except where noted.
 */
#ifndef PSVPP_HEURISTIC_PARAMS_H
#define PSVPP_HEURISTIC_PARAMS_H

/** @brief All tunable parameters of the heuristic. */
typedef struct AlnsParams {
    /* ---- multi-start (Algorithm 1) ---- */
    int starts;           /**< \f$\alpha\f$: number of initial solutions (50) */
    int start_iterations; /**< \f$\beta\f$: ALNS iterations applied to each start (100) */
    int iterations;       /**< \f$\gamma\f$: total number of ALNS iterations (300,000) */

    /* ---- ALNS (Algorithm 2) ---- */
    int ls_interval; /**< \f$\delta\f$: local search every this many iterations (50) */
    int segment;     /**< \f$\rho\f$: segment length for weight updates (50) */
    double t0;       /**< \f$\theta'_0\f$: initial temperature relative to the initial cost (1/50) */
    double tf;       /**< \f$\theta'_F\f$: final temperature relative to the initial cost (1/500) */

    /* ---- removal operators (Section 4.4) ---- */
    double removal_min; /**< \f$n^-_1\f$: lower bound of the removed fraction (0.10) */
    double removal_max; /**< \f$n^-_2\f$: upper bound of the removed fraction (0.25) */
    int removal_extra;  /**< installations removed on top of the drawn fraction (3, as in the original code) */
    double determinism; /**< \f$p\f$: level of determinism of the guided removals (4) */

    /* ---- insertion operators (Section 4.4) ---- */
    int insertion_sample; /**< removed installations examined per step by IO1 and IO3 (5) */

    /* ---- adaptive layer (Section 4.6) ---- */
    bool adaptive;       /**< roulette-wheel operator selection; false = uniform selection */
    double reaction;     /**< \f$\tau_A\f$: reaction factor (0.05) */
    double score_best;   /**< \f$\sigma_1\f$: score for a new best solution (5) */
    double score_better; /**< \f$\sigma_2\f$: score for improving the current solution (0) */
    double score_accept; /**< \f$\sigma_3\f$: score for an accepted worse solution (10) */

    /* ---- feasibility penalties (Section 4.2) ---- */
    double penalty_start;     /**< \f$\omega_0\f$: initial penalty weight (700) */
    double penalty_increase;  /**< factor applied when too many solutions are infeasible (1.15) */
    double penalty_decrease;  /**< factor applied when too few solutions are infeasible (0.90) */
    double target_infeasible; /**< target fraction of infeasible solutions (0.60) */
    double penalty_min;       /**< weights are not decreased below this value (10) */
    double penalty_max_load;  /**< capacity weight is not increased above this value (2500) */
    double penalty_max_time;  /**< time weight is not increased above this value (1500) */
    double penalty_max_berth; /**< departure weight is not increased above this value (2,000,000) */

    /* ---- voyages (Section 1.1, "Route duration") ---- */
    int max_route_size; /**< installations per voyage, at most (8); more count as a capacity violation */

    /* ---- local search (Algorithm 3) ---- */
    int ls_max_rounds; /**< maximum number of improvement steps per local search call (80) */
} AlnsParams;

/** @brief Returns the default parameter set (see the field documentation). */
AlnsParams alns_default_params(void);

#endif /* PSVPP_HEURISTIC_PARAMS_H */
