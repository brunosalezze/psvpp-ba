// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file local_search.c
 * @brief Voyage and pattern improvement neighbourhoods (see local_search.h).
 */
#include "heuristic/local_search.h"

#include "common/rng.h"
#include "heuristic/insertion.h"
#include "heuristic/moves.h"

/** Longest sequence moved by the exchange neighbourhood (paper: k in {0..3}). */
enum { MAX_EXCHANGE = 3 };

/** A sequence of consecutive visits detached from its voyage. */
typedef struct Sequence {
    int node[MAX_EXCHANGE]; /**< visits, in route order */
    int count;              /**< number of visits */
    int before;             /**< node preceding the sequence */
    int after;              /**< node following the sequence */
} Sequence;

/* State of one exchange attempt: the two sequences and the accumulated cost. */
typedef struct Exchange {
    Search *S;
    Solution *s;
    Sequence a, b;
    Cost delta;
} Exchange;

/* Inserts the visits of @p q between its before/after nodes (no-op when empty); adds the cost change to x. */
static void seq_insert(Exchange *x, Sequence *q)
{
    if (q->count > 0)
        x->delta += move_segment(x->S, x->s, q->node, q->count, q->before, q->after, MOVE_INSERT, true);
}

/* Removes the visits of @p q from between its before/after nodes; adds the cost change to x. */
static void seq_remove(Exchange *x, Sequence *q)
{
    if (q->count > 0)
        x->delta += move_segment(x->S, x->s, q->node, q->count, q->before, q->after, MOVE_REMOVE, true);
}

/* Reverses the order of the visits of @p q. */
static void seq_reverse(Sequence *q)
{
    for (int i = 0, j = q->count - 1; i < j; i++, j--) {
        const int t = q->node[i];
        q->node[i] = q->node[j];
        q->node[j] = t;
    }
}

/* Detaches the (up to) @p n visits following @p before (stopping at the voyage end). */
static void seq_detach(Exchange *x, Sequence *q, int before, int n)
{
    const Node *nd = x->s->node;
    const Instance *in = x->S->inst;
    q->before = before;
    q->count = 0;
    int u = nd[before].next;
    for (; q->count < n && !inst_is_head(in, u); u = nd[u].next)
        q->node[q->count++] = u;
    q->after = u;
    seq_remove(x, q);
}

/* Swaps the contents (not the positions) of the two sequences. */
static void seq_swap_contents(Exchange *x)
{
    const int lim = x->a.count > x->b.count ? x->a.count : x->b.count;
    const int t = x->a.count;
    x->a.count = x->b.count;
    x->b.count = t;
    for (int i = lim - 1; i >= 0; --i) {
        const int v = x->a.node[i];
        x->a.node[i] = x->b.node[i];
        x->b.node[i] = v;
    }
}

/* True if the (up to) n visits following @p before are all the visits of its voyage. */
static bool takes_whole_voyage(const Instance *in, const Solution *s, int before, int n)
{
    if (!inst_is_head(in, before))
        return false;
    int u = s->node[before].next, count = 0;
    while (count < n && !inst_is_head(in, u)) {
        u = s->node[u].next;
        count++;
    }
    return count > 0 && inst_is_head(in, u);
}

/* True if the two keys identify the same voyage. */
static bool same_voyage(VoyageKey a, VoyageKey b)
{
    return a.period == b.period && a.berth == b.berth && a.type == b.type && a.slot == b.slot;
}

/**
 * Exchanges the @p na visits after @p node_a with the @p nb visits after
 * @p node_b, trying the reversed sequences as well. Keeps the first variant
 * that improves the cost by more than the tolerance and returns true;
 * otherwise restores the solution and returns false.
 *
 * Restoring by the reverse moves is exact unless a voyage becomes empty on
 * the way: an emptied voyage rejoins its berth's departure queue at the end,
 * which changes the berth schedule. When a voyage can be emptied (a sequence
 * is the whole voyage, or both sequences are in one voyage that they may
 * cover), the solution is restored from a snapshot instead.
 */
static bool exchange_if_improving(Search *S, Solution *s, int node_a, int na, int node_b, int nb)
{
    Exchange x = {.S = S, .s = s, .delta = 0};
    const Cost tol = -(Cost)LS_IMPROVEMENT_TOLERANCE;
    const Instance *in = S->inst;
    const VoyageKey va = s->node[node_a].at, vb = s->node[node_b].at;
    const bool snapshot = same_voyage(va, vb)
                              ? sol_voyage_of(s, va)->size <= na + nb
                              : takes_whole_voyage(in, s, node_a, na) || takes_whole_voyage(in, s, node_b, nb);
    if (snapshot)
        sol_copy(S->undo, s);
    seq_detach(&x, &x.a, node_a, na);
    seq_detach(&x, &x.b, node_b, nb);
    seq_swap_contents(&x);
    seq_insert(&x, &x.a);
    seq_insert(&x, &x.b);
    if (x.delta < tol)
        return true;
    if (x.a.count > 1) {
        seq_remove(&x, &x.a);
        seq_reverse(&x.a);
        seq_insert(&x, &x.a);
        if (x.delta < tol)
            return true;
    }
    if (x.b.count > 1) {
        seq_remove(&x, &x.b);
        seq_reverse(&x.b);
        seq_insert(&x, &x.b);
        if (x.delta < tol)
            return true;
    }
    if (x.a.count > 1 && x.b.count > 1) { /* both sequences reversed back to their original orientation */
        seq_remove(&x, &x.a);
        seq_reverse(&x.a);
        seq_insert(&x, &x.a);
        if (x.delta < tol)
            return true;
    }
    if (snapshot) {
        sol_copy(s, S->undo);
        return false;
    }
    /* Undo: take both sequences out, restore their orientation and put them back. */
    seq_remove(&x, &x.a);
    seq_remove(&x, &x.b);
    if (x.a.count > 1 && x.b.count > 1)
        seq_reverse(&x.b); /* a is already in its original orientation */
    else {
        if (x.a.count > 1)
            seq_reverse(&x.a);
        if (x.b.count > 1)
            seq_reverse(&x.b);
    }
    seq_swap_contents(&x);
    seq_insert(&x, &x.a);
    seq_insert(&x, &x.b);
    return false;
}

/**
 * Tries the exchanges between visit @p u (with its predecessor a, successor x)
 * and position @p v (with predecessor b, successor y), in random order.
 * Exchange m (0..7) moves (cnt_a visits after a) <-> (cnt_b visits after anchor):
 *   0:(1,v,0) 1:(2,v,0) 2:(3,v,0) 3:(1,b,1) 4:(2,b,1) 5:(2,b,2) 6:(3,b,1) 7:(3,b,2)
 * Combinations that would overlap or run past the voyage are skipped.
 */
static bool try_exchanges(Search *S, Solution *s, int u, int v)
{
    const Instance *in = S->inst;
    const Node *nd = s->node;
    static const struct {
        int count_a, anchor_is_b, count_b;
    } MOVES[8] = {
        {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {1, 1, 1}, {2, 1, 1}, {2, 1, 2}, {3, 1, 1}, {3, 1, 2},
    };
    int order[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    int num = 8;
    const int y = nd[v].next, x = nd[u].next, b = nd[v].prev, a = nd[u].prev;

    if (u == y)
        return false;
    if (inst_is_head(in, v) || x == v) {
        num = 1;
    } else if (inst_is_head(in, x) || x == b) {
        num = 4;
    } else {
        const int k = nd[x].next;
        if (inst_is_head(in, y) || y == a) {
            num = 7;
            order[5] = -1;
        }
        if (inst_is_head(in, k) || k == b)
            num = 6;
    }
    if (!inst_is_head(in, x) && x != v) {
        const int k = nd[x].next;
        if (inst_is_head(in, k) || k == v)
            order[2] = -1;
    }
    rng_shuffle_int(&S->rng, order, num);
    for (int i = 0; i < num; i++) {
        if (order[i] < 0)
            continue;
        const int m = order[i];
        if (exchange_if_improving(S, s, a, MOVES[m].count_a, MOVES[m].anchor_is_b ? b : v, MOVES[m].count_b))
            return true;
    }
    return false;
}

bool ls_voyage_improvement(Search *S, Solution *s)
{
    const Instance *in = S->inst;
    int *order = S->ls_order, *cand = S->ls_cand;
    for (int i = 0; i < in->num_visits; i++)
        order[i] = i;
    rng_shuffle_int(&S->rng, order, in->num_visits);
    for (int i = 0; i < in->num_visits; i++) {
        const int u = order[i];
        const int n = insertion_candidates(s, in->visit_inst[u], s->node[u].at.period, cand);
        rng_shuffle_int(&S->rng, cand, n);
        for (int j = 0; j < n; j++)
            if (try_exchanges(S, s, u, cand[j]))
                return true;
    }
    return false;
}

bool ls_pattern_improvement(Search *S, Solution *s)
{
    const Instance *in = S->inst;
    int *order = S->ls_order;
    for (int i = 0; i < in->num_installations; i++)
        order[i] = i;
    rng_shuffle_int(&S->rng, order, in->num_installations);

    PatternChoice best = {.before = S->choice_buf};
    for (int i = 0; i < in->num_installations; i++) {
        const int inst = order[i];
        best.inst = inst;
        sol_copy(S->undo, s); /* restored by copy, see exchange_if_improving() */
        const Cost saving = -insertion_remove_installation(S, s, inst);
        insertion_best_pattern(S, s, &best, 0);
        const Cost cost = insertion_apply(S, s, &best);
        if (saving - cost > LS_IMPROVEMENT_TOLERANCE)
            return true;
        sol_copy(s, S->undo);
    }
    return false;
}

void local_search(Search *S, Solution *s)
{
    bool improved[2] = {true, true};
    int which = rng_int(&S->rng, 0, 1);
    for (int round = 0; round < S->par->ls_max_rounds && (improved[0] || improved[1]); round++) {
        improved[which] = which ? ls_voyage_improvement(S, s) : ls_pattern_improvement(S, s);
        which = 1 - which;
    }
    search_update_cost(S, s);
}
