// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file test_moves.c
 * @brief Unit tests of the incremental evaluation of the heuristic moves
 *        (`src/heuristic/moves.c`).
 *
 * Two properties of the evaluator are checked here. Both depend on the exact
 * floating-point model of the build (`-ffp-contract=off -fno-fast-math` in the
 * Makefile) and neither survives a change of that model:
 *
 *  1. an insertion evaluated with a pruning bound that its cost does not reach
 *     must give exactly the same ::Cost as the same insertion evaluated without
 *     a bound: the insertion operator (`insertion.c`) compares these sums of
 *     floats with each other to choose where a visit goes, so a bound may only
 *     skip an evaluation, never change the value of the ones it completes;
 *  2. a removal followed by the same insertion cancels exactly: the two cost
 *     deltas add to zero and the solution comes back bit for bit, which is how
 *     the local search undoes a move that did not improve (see
 *     `exchange_if_improving()`), and the delta of a move agrees with what a
 *     full recomputation of the penalised objective gives.
 *
 * usage: build/tests/test_moves [instances-directory]
 */
#include "common/instance.h"
#include "common/util.h"
#include "heuristic/alns.h"
#include "heuristic/construction.h"
#include "heuristic/insertion.h"
#include "heuristic/moves.h"
#include "heuristic/search.h"
#include "heuristic/solution.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** The smallest instance of the set, and a larger one with two berths. */
static const char *const tested_instances[] = {"p10-1-1", "p30-2-1"};

/** A fixed seed makes every constructed solution, and every check, reproducible. */
static const uint64_t test_seed = 1;

static int checks;   /**< checks run so far */
static int failures; /**< checks that failed */

/** Reports a failed check with the expression that did not hold. */
#define CHECK(condition)                                                                  \
    do {                                                                                  \
        checks++;                                                                         \
        if (!(condition)) {                                                               \
            failures++;                                                                   \
            fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
        }                                                                                 \
    } while (0)

/** Fails unless the two costs are the same float, bit for bit (0.0 and -0.0 differ). */
#define CHECK_SAME(actual, expected)                                                               \
    do {                                                                                           \
        checks++;                                                                                  \
        const Cost actual_ = (actual), expected_ = (expected);                                     \
        if (!cost_equal(actual_, expected_)) {                                                     \
            failures++;                                                                            \
            fprintf(stderr, "%s:%d: %s is %.9g and %s is %.9g, but they must be the same float\n", \
                    __FILE__, __LINE__, #actual, (double)actual_, #expected, (double)expected_);   \
        }                                                                                          \
    } while (0)

/** Fails unless the two costs are at most @p tolerance apart. */
#define CHECK_CLOSE(actual, expected, tolerance)                                                      \
    do {                                                                                              \
        checks++;                                                                                     \
        const double actual_ = (double)(actual), expected_ = (double)(expected);                      \
        const double diff_ = actual_ > expected_ ? actual_ - expected_ : expected_ - actual_;         \
        if (!(diff_ <= (double)(tolerance))) {                                                        \
            failures++;                                                                               \
            fprintf(stderr, "%s:%d: %s is %.9g and %s is %.9g, more than %.6g apart\n",               \
                    __FILE__, __LINE__, #actual, actual_, #expected, expected_, (double)(tolerance)); \
        }                                                                                             \
    } while (0)

/** True if the two costs are the same float, bit for bit. */
static bool cost_equal(Cost a, Cost b) { return memcmp(&a, &b, sizeof(Cost)) == 0; }

/**
 * @brief Tolerance of a cost comparison that cannot be exact.
 *
 * The penalised objective is a sum of floats; the incremental delta of a move and a full
 * recomputation add the same terms in a different order, so they may differ by a rounding.
 */
static double cost_tolerance(Cost reference)
{
    const double size = fabs((double)reference);
    return 1e-4 * (size > 1.0 ? size : 1.0);
}

/**
 * @brief Looks an instance file up in the usual places.
 *
 * @param dir  directory given on the command line, or NULL to search the default ones.
 * @param name instance name, without a directory.
 * @param out  receives the path.
 * @param size size of @p out.
 * @return @p out when the file was found, NULL otherwise.
 */
static char *find_instance(const char *dir, const char *name, char *out, size_t size)
{
    static const char *const places[] = {"instances", "../instances", "../../instances"};
    for (unsigned i = 0; i <= sizeof(places) / sizeof(*places); i++) {
        const char *base = i == 0 ? dir : places[i - 1];
        if (base == NULL)
            continue;
        snprintf(out, size, "%s/%s", base, name);
        FILE *const f = fopen(out, "r");
        if (f != NULL) {
            fclose(f);
            return out;
        }
    }
    return NULL;
}

/** Number of visits the voyages of @p s serve. */
static int routed_visits(const Solution *s)
{
    const Instance *in = s->inst;
    int n = 0;
    for (int head = in->num_visits; head < in->num_nodes; head++)
        for (int u = s->node[head].next; u != head; u = s->node[u].next)
            n++;
    return n;
}

/**
 * @brief Collects the visits that can be taken out of their voyage without emptying it.
 *
 * A voyage that becomes empty leaves the queue of its berth day and comes back at its end,
 * which changes the berth schedule: only a visit of a voyage with two visits or more is
 * restored exactly by the reverse move.
 *
 * @param s      solution.
 * @param visits receives one node per pair (at least ::Instance::num_visits entries).
 * @param before receives the node preceding each of them (same length as @p visits).
 * @return Number of pairs written.
 */
static int movable_visits(const Solution *s, int *visits, int *before)
{
    const Instance *in = s->inst;
    int n = 0;
    for (int head = in->num_visits; head < in->num_nodes; head++) {
        int size = 0;
        for (int u = s->node[head].next; u != head; u = s->node[u].next)
            size++;
        if (size < 2)
            continue;
        for (int u = s->node[head].next; u != head; u = s->node[u].next) {
            visits[n] = u;
            before[n] = s->node[u].prev;
            n++;
        }
    }
    return n;
}

/**
 * @brief A removal and the insertion of the same visit at the same place cancel exactly.
 *
 * Checks the cost deltas, the penalised cost recomputed from scratch, and every byte of the
 * solution block: the linked lists, the voyages, the berth days and the fleet counters.
 */
static void test_removal_round_trip(Search *S, Solution *s)
{
    const Instance *in = S->inst;
    int *visits = XMALLOC_ARRAY(int, in->num_visits);
    int *before = XMALLOC_ARRAY(int, in->num_visits);
    Solution *snapshot = sol_create(in);
    search_update_cost(S, s);

    const int n = movable_visits(s, visits, before);
    CHECK(n > 0);
    for (int i = 0; i < n; i++) {
        const int u = visits[i];
        sol_copy(snapshot, s);
        const Cost removed = move_remove(S, s, u);
        const Cost inserted = move_insert(S, s, u, before[i]);

        CHECK_SAME(removed + inserted, 0.0f);
        CHECK(memcmp(s->block, snapshot->block, s->block_size) == 0);
        search_update_cost(S, s);
        CHECK_SAME(s->cost, snapshot->cost);
        CHECK(s->infeasible == snapshot->infeasible);
        CHECK(sol_objective(s) == sol_objective(snapshot));

        /* the removal on its own must move the recomputed cost by the delta it announced */
        move_remove(S, s, u);
        search_update_cost(S, s);
        CHECK_CLOSE(s->cost, snapshot->cost + removed, cost_tolerance(snapshot->cost));
        sol_copy(s, snapshot);
    }
    sol_free(snapshot);
    free(before);
    free(visits);
}

/**
 * @brief A pruning bound that is not reached never changes the cost an insertion reports.
 *
 * Every installation is taken out of the solution and each position the insertion operator
 * would look at is evaluated with and without a bound; then the installation goes back.
 */
static void test_insertion_bounds(Search *S, Solution *s)
{
    const Instance *in = S->inst;
    int *candidates = XMALLOC_ARRAY(int, in->num_installations + in->num_vessels);
    PatternChoice choice = {.inst = 0, .before = S->choice_buf};
    search_update_cost(S, s);
    int evaluated = 0, pruned = 0, negative = 0;

    for (int inst = 0; inst < in->num_installations; inst++) {
        const Installation *u = &in->inst[inst];
        /* the visits of one installation become the unrouted visits to insert back */
        insertion_remove_installation(S, s, inst);
        CHECK(routed_visits(s) == in->num_visits - u->frequency);
        for (int p = 0; p < in->num_periods; p++) {
            const int node = u->first_visit + p % u->frequency;
            const int n = insertion_candidates(s, inst, p, candidates);
            CHECK(n > 0);
            for (int i = 0; i < n; i++) {
                const int before = candidates[i];
                const Cost exact = move_insert_cost(S, s, node, before);
                evaluated++;
                negative += exact < 0.0f;

                /* no bound: the same evaluation through the bounded entry point */
                CHECK_SAME(move_insert_cost_below(S, s, node, before, INFINITY), exact);
                /* a bound no cost can reach: the exact cost again, with the same bits */
                const Cost far = exact + fabsf(exact) + (Cost)COST_SCALE;
                CHECK_SAME(move_insert_cost_below(S, s, node, before, far), exact);
                /* a bound below the cost: the evaluation may stop, but never below the bound */
                const Cost bound = exact * 0.5f;
                const Cost stopped = move_insert_cost_below(S, s, node, before, bound);
                CHECK(stopped >= bound || cost_equal(stopped, exact));
                CHECK(stopped <= exact + cost_tolerance(exact));
                pruned += !cost_equal(stopped, exact);
            }
        }
        /* put the installation back, as the insertion operator would, before the next one */
        choice.inst = inst;
        insertion_best_pattern(S, s, &choice, 0);
        insertion_apply(S, s, &choice);
        CHECK(routed_visits(s) == in->num_visits);
        CHECK(choice.pattern >= 0 && choice.pattern < u->num_patterns);
        CHECK(s->pattern[inst] == choice.pattern);
    }
    search_update_cost(S, s);
    free(candidates);
    printf("   %d insertion positions evaluated: %d stopped at the bound, %d cheaper than doing nothing\n",
           evaluated, pruned, negative);
}

/** Runs the checks of one instance: on a random start, then on a solution the search improved. */
static void test_instance(const char *path)
{
    Instance *in = inst_load(path);
    const AlnsParams par = alns_default_params(); /* the tuned parameters of the paper */
    Search S;
    search_init(&S, in, &par, test_seed);
    Solution *s = sol_create(in);

    printf("  %s: %d installations, %d weekly visits\n", in->name, in->num_installations, in->num_visits);
    construct_random_solution(&S, s);
    CHECK(routed_visits(s) == in->num_visits);
    test_insertion_bounds(&S, s);
    test_removal_round_trip(&S, s);

    /* the same checks on a solution the search has actually worked on */
    AlnsReport report = {0}; /* alns_run() accumulates into the report it is given */
    sol_clear(s);
    construct_random_solution(&S, s);
    search_update_cost(&S, s);
    alns_run(&S, s, 500, &report);
    CHECK(report.iterations == 500);
    CHECK(routed_visits(s) == in->num_visits);
    test_removal_round_trip(&S, s);
    test_insertion_bounds(&S, s);

    sol_free(s);
    search_free(&S);
    inst_free(in);
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : getenv("PSVPP_INSTANCES");
    char path[4096];
    for (unsigned i = 0; i < sizeof(tested_instances) / sizeof(*tested_instances); i++) {
        if (find_instance(dir, tested_instances[i], path, sizeof path) == NULL) {
            fprintf(stderr, "%s: instance not found (give the instances directory as the first argument)\n",
                    tested_instances[i]);
            return EXIT_FAILURE;
        }
        test_instance(path);
    }
    if (failures > 0) {
        fprintf(stderr, "test_moves: %d of %d checks failed\n", failures, checks);
        return EXIT_FAILURE;
    }
    printf("test_moves: %d checks passed\n", checks);
    return EXIT_SUCCESS;
}
