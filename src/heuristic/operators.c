// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file operators.c
 * @brief Removal and insertion operators (see operators.h).
 */
#include "heuristic/operators.h"

#include "common/util.h"
#include "heuristic/insertion.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/** "No candidate yet" value of the sampled insertions (larger than any insertion cost). */
static const Cost NO_CHOICE = 400000000000.0f;
/** Offset guarding the ranking index against floating-point truncation. */
static const float INDEX_ROUNDING = 0.1f;

void pool_init(InstallationPool *pool, const Instance *inst)
{
    const int N = inst->num_installations;
    pool->removed = XMALLOC_ARRAY(int, N);
    pool->routed = XMALLOC_ARRAY(int, N);
    pool->key = XCALLOC_ARRAY(int, N);
    for (int i = 0; i < N; i++)
        pool->routed[i] = i;
    pool->num_routed = N;
    pool->num_removed = 0;
}

void pool_free(InstallationPool *pool)
{
    free(pool->removed);
    free(pool->routed);
    free(pool->key);
    *pool = (InstallationPool){0};
}

/* Moves routed[idx] to the removed list, keeping the order of the routed list. */
static void pool_take(InstallationPool *pool, int idx)
{
    pool->removed[pool->num_removed++] = pool->routed[idx];
    pool->num_routed--;
    memmove(&pool->routed[idx], &pool->routed[idx + 1], sizeof(int) * (size_t)(pool->num_routed - idx));
}

/* Moves removed[idx] back to the end of the routed list. */
static void pool_put_back(InstallationPool *pool, int idx)
{
    pool->routed[pool->num_routed++] = pool->removed[idx];
    pool->removed[idx] = pool->removed[--pool->num_removed];
}

/* ------------------------------------------------------------------------- */
/*  Removal                                                                  */
/* ------------------------------------------------------------------------- */

/* Sorts the routed installations by decreasing ranking key. The sort is stable (installations
 * with equal keys keep their order), so a run does not depend on the C library's qsort; an
 * insertion sort is enough for lists of this size, which are also partly sorted already. */
static void sort_routed_by_key(InstallationPool *pool)
{
    int *a = pool->routed;
    const int *key = pool->key;
    for (int i = 1; i < pool->num_routed; i++) {
        const int x = a[i], kx = key[x];
        int j = i;
        for (; j > 0 && key[a[j - 1]] < kx; j--)
            a[j] = a[j - 1];
        a[j] = x;
    }
}

/* Shaw relatedness (paper, Eq. 37) of every routed installation to @p seed, as an integer key. */
static void shaw_keys(const Instance *in, InstallationPool *pool, int seed, float w_dist, float w_demand, float w_freq)
{
    const Installation *r = &in->inst[seed];
    const float max_dist = (float)in->max_dist, max_demand = (float)in->max_demand;
    for (int i = 0; i < pool->num_routed; i++) {
        const int j = pool->routed[i];
        const Installation *u = &in->inst[j];
        const float p1 = w_dist * (float)inst_dist(in, seed, j) / max_dist;
        const float p2 = w_demand * (float)iabs(r->demand - u->demand) / max_demand;
        const float p3 = w_freq * (float)iabs(u->frequency - r->frequency) / (float)in->num_periods;
        pool->key[j] = (int)(1e6f * (p1 + p2 + p3));
    }
}

void apply_removal(Search *S, Solution *s, InstallationPool *pool, RemovalOp op)
{
    const Instance *in = S->inst;
    const AlnsParams *par = S->par;
    const bool guided = op != RO4_RANDOM;

    switch (op) {
    case RO1_SHAW:
        shaw_keys(in, pool, pool->routed[rng_int(&S->rng, 0, pool->num_routed - 1)], 1.0f, 0.4f, 0.4f);
        break;
    case RO2_SHAW_DISTANCE: {
        const int seed = pool->routed[rng_int(&S->rng, 0, pool->num_routed - 1)];
        for (int i = 0; i < pool->num_routed; i++)
            pool->key[pool->routed[i]] = inst_dist(in, seed, pool->routed[i]);
        break;
    }
    case RO3_SHAW_FREQUENCY:
        for (int i = 0; i < pool->num_routed; i++) {
            const int j = pool->routed[i];
            const int f = in->inst[j].frequency;
            pool->key[j] = (int)(1000.0 * fabs((float)in->num_periods / 2 - (float)f));
        }
        break;
    case RO5_WORST:
        for (int i = 0; i < pool->num_routed; i++) {
            const int j = pool->routed[i];
            pool->key[j] = (int)insertion_removal_cost(S, s, j);
        }
        break;
    case RO4_RANDOM:
    case NUM_REMOVAL_OPS: break;
    }
    if (guided)
        sort_routed_by_key(pool);

    const int n = pool->num_routed;
    int count = rng_int(&S->rng, (int)(par->removal_min * n), (int)(par->removal_max * n)) + par->removal_extra;
    count = imin(count, n);
    for (; count > 0; count--) {
        int idx;
        if (guided) {
            const double a = pow(rng_unit(&S->rng), par->determinism);
            idx = (int)(INDEX_ROUNDING + (1 - a) * (pool->num_routed - 1));
        } else {
            idx = rng_int(&S->rng, 0, pool->num_routed - 1);
        }
        insertion_remove_installation(S, s, pool->routed[idx]);
        pool_take(pool, idx);
    }
}

/* ------------------------------------------------------------------------- */
/*  Insertion                                                                */
/* ------------------------------------------------------------------------- */

/* Inserts, among (up to) AlnsParams::insertion_sample removed installations
 * taken in random order, the one with the smallest value (cost or regret). */
static void insert_best_of_sample(Search *S, Solution *s, InstallationPool *pool, int regret_k)
{
    const int P = S->inst->num_periods;
    PatternChoice best = {.before = S->choice_buf}, cand = {.before = S->choice_buf + P};
    do {
        best.value = NO_CHOICE;
        int best_idx = 0;
        rng_shuffle_int(&S->rng, pool->removed, pool->num_removed);
        const int limit = imin(S->par->insertion_sample, pool->num_removed);
        for (int i = 0; i < limit; i++) {
            cand.inst = pool->removed[i];
            insertion_best_pattern(S, s, &cand, regret_k);
            if (cand.value < best.value || i == 0) {
                memcpy(best.before, cand.before, sizeof(int) * (size_t)P);
                best.inst = cand.inst;
                best.pattern = cand.pattern;
                best.value = cand.value;
                best_idx = i;
            }
        }
        insertion_apply(S, s, &best);
        pool_put_back(pool, best_idx);
    } while (pool->num_removed > 0);
}

void apply_insertion(Search *S, Solution *s, InstallationPool *pool, InsertionOp op)
{
    switch (op) {
    case IO1_DEEP_GREEDY: insert_best_of_sample(S, s, pool, 0); break;
    case IO3_REGRET: insert_best_of_sample(S, s, pool, rng_int(&S->rng, 2, 4)); break;
    case IO2_GREEDY: {
        PatternChoice c = {.before = S->choice_buf};
        do {
            const int idx = rng_int(&S->rng, 0, pool->num_removed - 1);
            c.inst = pool->removed[idx];
            insertion_best_pattern(S, s, &c, 0);
            insertion_apply(S, s, &c);
            pool_put_back(pool, idx);
        } while (pool->num_removed > 0);
        break;
    }
    case NUM_INSERTION_OPS: break;
    }
}

const char *removal_op_name(RemovalOp op)
{
    static const char *const names[NUM_REMOVAL_OPS] = {"RO1-shaw", "RO2-shaw-distance", "RO3-shaw-frequency",
                                                       "RO4-random", "RO5-worst"};
    return op < NUM_REMOVAL_OPS ? names[op] : "?";
}

const char *insertion_op_name(InsertionOp op)
{
    static const char *const names[NUM_INSERTION_OPS] = {"IO1-deep-greedy", "IO2-greedy", "IO3-regret"};
    return op < NUM_INSERTION_OPS ? names[op] : "?";
}
