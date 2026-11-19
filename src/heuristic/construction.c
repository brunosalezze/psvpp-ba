// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file construction.c
 * @brief Fix and Split construction (see construction.h).
 */
#include "heuristic/construction.h"

#include "common/util.h"
#include "heuristic/moves.h"

#include <stdlib.h>
#include <string.h>

static const float SPLIT_START_FILL = 0.75f; /**< initial fill fraction of the vessels */
static const float SPLIT_FILL_STEP = 0.25f;  /**< fill increase after a failed split */
enum { SPLIT_ATTEMPTS = 10 };                /**< fill levels tried before giving up */

/* Fix phase: random pattern per installation (written to @p pattern) and shuffled giant tour per
 * day (Search::tour, Search::tour_size). */
static void fix_patterns(Search *S, int8_t *pattern)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, N = in->num_installations;
    for (int p = 0; p < P; p++)
        S->tour_size[p] = 0;
    for (int i = 0; i < N; i++) {
        const Installation *u = &in->inst[i];
        const int pat = rng_int(&S->rng, 0, u->num_patterns - 1);
        pattern[i] = (int8_t)pat;
        for (int p = 0, k = 0; p < P; p++)
            if (inst_pattern_visits(u, P, pat, p))
                S->tour[p * N + S->tour_size[p]++] = (u->first_visit + k++);
    }
    for (int p = 0; p < P; p++)
        rng_shuffle_int(&S->rng, &S->tour[p * N], S->tour_size[p]);
}

/* Split of one day's giant tour. Returns true if every visit was routed. */
static bool split_day(Search *S, Solution *s, int p, float fill, const int *tour, int n)
{
    const Instance *in = S->inst;
    int free_slots = 0;
    for (int b = 0; b < in->num_berths; b++)
        for (int k = 0; k < in->num_types; k++)
            for (int l = 0; l < in->type[k].max_slots; l++)
                free_slots += sol_voyage(s, inst_vessel(in, b, k, l), p)->size == 0;

    int i = 0;
    while (free_slots > 0) {
        const int k = rng_int(&S->rng, 0, in->num_types - 1);
        const int b = rng_int(&S->rng, 0, in->num_berths - 1);
        int head = -1;
        for (int l = 0; l < in->type[k].max_slots; l++) {
            const int v = inst_vessel(in, b, k, l);
            if (sol_voyage(s, v, p)->size == 0) {
                head = inst_head_node(in, p, v);
                break;
            }
        }
        if (head >= 0) {
            const int limit = (int)((float)in->type[k].capacity * fill);
            int j = i, load = 0;
            for (; load < limit && j < n && j - i < S->par->max_route_size; j++)
                load += in->visit_load[tour[j]];
            if (j > i) {
                move_insert_many(S, s, &tour[i], j - i, head);
                i = j;
                free_slots--;
            }
        }
        if (i == n)
            return true;
    }
    return false;
}

void construct_random_solution(Search *S, Solution *s)
{
    const Instance *in = S->inst;
    const int P = in->num_periods, N = in->num_installations;
    int8_t *pattern = XMALLOC_ARRAY(int8_t, N); /* sol_clear() wipes the patterns of the solution */
    fix_patterns(S, pattern);

    float fill = SPLIT_START_FILL;
    for (int attempt = 0; attempt < SPLIT_ATTEMPTS; attempt++, fill += SPLIT_FILL_STEP) {
        sol_clear(s);
        int failures = 0;
        for (int p = 0; p < P; p++)
            if (S->tour_size[p] > 0)
                failures += !split_day(S, s, p, fill, &S->tour[p * N], S->tour_size[p]);
        if (failures == 0) {
            sol_compact(s);
            sol_check_feasibility(s);
            search_update_cost(S, s);
            memcpy(s->pattern, pattern, sizeof(int8_t) * (size_t)N);
            free(pattern);
            return;
        }
    }
    die("%s: the vessel slots cannot carry the weekly demand of some day; "
        "increase the slots per berth of the vessel types in the instance file",
        in->name);
}
