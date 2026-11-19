// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file insertion.c
 * @brief Best insertion over delivery patterns (see insertion.h).
 */
#include "heuristic/insertion.h"

#include "heuristic/moves.h"

#include <math.h>

/** Initial value of "best so far" comparisons (larger than any realistic insertion cost). */
static const Cost NO_COST = 3999999999.0f;

int insertion_candidates(const Solution *s, int inst, int period, int *out)
{
    const Instance *in = s->inst;
    const int N = in->num_installations, limit = in->num_neighbors[inst];
    const int *nb = &in->neighbors[inst * (N - 1)];
    const int8_t *served = &s->visit_count[period * N];
    const int *node = &s->visit_node[period * N];
    int n = 0;
    /* Visits of the nearest neighbours served on this day. The slot is always written and kept
     * only when the neighbour is served: this unpredictable test costs more as a branch. */
    for (int i = 0; i < N - 1 && n < limit; i++) {
        out[n] = node[nb[i]];
        n += served[nb[i]] != 0;
    }
    /* Heads of the active vessels, and of the first inactive slot of each (berth, type). */
    for (int b = in->num_berths - 1; b >= 0; --b)
        for (int k = in->num_types - 1; k >= 0; --k) {
            const int active = *sol_active(s, b, k);
            for (int l = active - 1; l >= 0; --l)
                out[n++] = inst_head_node(in, period, inst_vessel(in, b, k, l));
            if (active < in->type[k].max_slots)
                out[n++] = inst_head_node(in, period, inst_vessel(in, b, k, active));
        }
    return n;
}

/** Largest regret level of the k-regret insertion (IO3 draws k in 2..4). */
enum { MAX_REGRET = 4 };

/* Best position of one visit of @p inst on @p period; with regret_k >= 2 the
 * returned cost is the regret of the day instead of the best cost.
 *
 * Only the regret_k cheapest costs matter (one for the greedy insertion), so the largest of
 * them bounds the evaluation of the remaining candidates, see move_insert_cost_below(). */
static InsertOption best_position(Search *S, Solution *s, int inst, int period, int regret_k)
{
    const int visit = S->inst->inst[inst].first_visit;
    const int n = insertion_candidates(s, inst, period, S->cand);
    const int k = regret_k > 1 ? (regret_k < MAX_REGRET ? regret_k : MAX_REGRET) : 1;
    InsertOption best = {.before = n > 0 ? S->cand[0] : -1, .cost = NO_COST};
    Cost low[MAX_REGRET] = {0}; /* the (up to) k cheapest costs found, in increasing order */
    int m = 0;
    for (int i = 0; i < n; i++) {
        const Cost bound = k == 1 ? best.cost : m < k ? INFINITY : low[k - 1];
        const Cost cost = move_insert_cost_below(S, s, visit, S->cand[i], bound);
        if (best.cost > cost)
            best = (InsertOption){.before = S->cand[i], .cost = cost};
        if (k > 1 && (m < k || cost < low[k - 1])) {
            int j = m < k ? m++ : k - 1;
            for (; j > 0 && low[j - 1] > cost; j--)
                low[j] = low[j - 1];
            low[j] = cost;
        }
    }
    if (k > 1) {
        Cost regret = 0;
        for (int i = 1; i < m; i++)
            regret += best.cost - low[i];
        best.cost = regret;
    }
    return best;
}

void insertion_best_pattern(Search *S, Solution *s, PatternChoice *choice, int regret_k)
{
    const Instance *in = S->inst;
    const int P = in->num_periods;
    const Installation *u = &in->inst[choice->inst];
    InsertOption *day = S->day_best;
    for (int p = 0; p < P; p++) {
        const InsertOption o = best_position(S, s, choice->inst, p, regret_k);
        day[p] = o.cost < NO_COST ? o : (InsertOption){.before = o.before, .cost = NO_COST};
        choice->before[p] = day[p].before;
    }
    Cost *pc = S->pattern_cost;
    for (int i = 0; i < u->num_patterns; i++) {
        pc[i] = 0;
        for (int p = 0; p < P; p++)
            if (inst_pattern_visits(u, P, i, p))
                pc[i] += day[p].cost;
    }
    choice->pattern = 0;
    for (int i = 1; i < u->num_patterns; i++)
        if (pc[choice->pattern] > pc[i])
            choice->pattern = i;
    choice->value = pc[choice->pattern];
}

Cost insertion_apply(Search *S, Solution *s, const PatternChoice *choice)
{
    const Instance *in = S->inst;
    const Installation *u = &in->inst[choice->inst];
    Cost total = 0;
    for (int p = 0, k = 0; p < in->num_periods; p++)
        if (inst_pattern_visits(u, in->num_periods, choice->pattern, p))
            total += move_insert(S, s, (u->first_visit + k++), choice->before[p]);
    s->pattern[choice->inst] = (int8_t)choice->pattern;
    return total;
}

Cost insertion_remove_installation(Search *S, Solution *s, int inst)
{
    const Installation *u = &S->inst->inst[inst];
    Cost total = 0;
    for (int k = 0; k < u->frequency; k++)
        total += move_remove(S, s, (u->first_visit + k));
    s->pattern[inst] = -1;
    return total;
}

Cost insertion_removal_cost(Search *S, Solution *s, int inst)
{
    const Installation *u = &S->inst->inst[inst];
    Cost total = 0;
    for (int k = 0; k < u->frequency; k++)
        total += move_remove_cost(S, s, (u->first_visit + k));
    return total;
}
