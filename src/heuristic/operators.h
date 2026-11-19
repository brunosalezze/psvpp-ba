// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file operators.h
 * @brief Removal (destroy) and insertion (repair) operators of the ALNS
 *        (paper, Section 4.4).
 *
 * A removal operator takes \f$n^- \in [n^-_1 |N^{OI}|, n^-_2 |N^{OI}|]\f$
 * installations (plus AlnsParams::removal_extra) out of the solution, all
 * their visits included. An insertion operator then reinserts every removed
 * installation with insertion_best_pattern().
 *
 * The guided removals rank the routed installations by a key and pick the
 * entry at position \f$\lfloor (1 - a^p)(n - 1) \rfloor\f$ of the ranking
 * sorted by decreasing key, with \f$a \in [0,1)\f$ uniform and \f$p\f$ the
 * level of determinism: low keys (very related or very costly installations)
 * are the most likely to be chosen.
 */
#ifndef PSVPP_HEURISTIC_OPERATORS_H
#define PSVPP_HEURISTIC_OPERATORS_H

#include "heuristic/search.h"
#include "heuristic/solution.h"

/** @brief Removal operators. */
typedef enum RemovalOp {
    RO1_SHAW,           /**< Shaw relatedness to a random seed installation */
    RO2_SHAW_DISTANCE,  /**< Shaw relatedness using distance only (nearest neighbours) */
    RO3_SHAW_FREQUENCY, /**< Shaw relatedness dominated by the visit frequency */
    RO4_RANDOM,         /**< uniformly random installations */
    RO5_WORST,          /**< installations whose removal saves the most */
    NUM_REMOVAL_OPS
} RemovalOp;

/** @brief Insertion operators. */
typedef enum InsertionOp {
    IO1_DEEP_GREEDY, /**< repeatedly insert the cheapest of a sample of removed installations */
    IO2_GREEDY,      /**< insert removed installations in random order at their best positions */
    IO3_REGRET,      /**< repeatedly insert the sampled installation with the largest k-regret */
    NUM_INSERTION_OPS
} InsertionOp;

/** @brief Installations currently removed from, and routed in, the solution. */
typedef struct InstallationPool {
    int *removed;    /**< [num_installations] removed installations \f$L^-\f$ */
    int num_removed; /**< size of @c removed */
    int *routed;     /**< [num_installations] installations in the solution */
    int num_routed;  /**< size of @c routed */
    int *key;        /**< [num_installations] ranking keys of the guided removals */
} InstallationPool;

/** @brief Initialises a pool in which every installation is routed. */
void pool_init(InstallationPool *pool, const Instance *inst);

/** @brief Releases the arrays of a pool. */
void pool_free(InstallationPool *pool);

/**
 * @brief Applies a removal operator to @p s.
 *
 * @param S    search context.
 * @param s    solution; on return the removed installations are unrouted.
 * @param pool on input every installation is routed; on output the removed
 *             ones are listed in InstallationPool::removed.
 * @param op   operator to apply.
 */
void apply_removal(Search *S, Solution *s, InstallationPool *pool, RemovalOp op);

/**
 * @brief Applies an insertion operator: reinserts every removed installation.
 *
 * @param S    search context.
 * @param s    solution.
 * @param pool on output every installation is routed again.
 * @param op   operator to apply.
 */
void apply_insertion(Search *S, Solution *s, InstallationPool *pool, InsertionOp op);

/** @brief Short name of a removal operator (for reports). */
const char *removal_op_name(RemovalOp op);

/** @brief Short name of an insertion operator (for reports). */
const char *insertion_op_name(InsertionOp op);

#endif /* PSVPP_HEURISTIC_OPERATORS_H */
