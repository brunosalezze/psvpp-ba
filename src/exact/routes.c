// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file routes.c
 * @brief Voyage enumeration with exact TSP sequencing (see routes.h).
 */
#include "exact/routes.h"

#include "common/timer.h"
#include "common/util.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Working state of the enumeration. */
typedef struct Builder {
    RoutePool *pool;
    const Instance *in;
    int depot;        /* location of the supply base (first berth) */
    int max_capacity; /* largest vessel capacity */
    int k;            /* installations in the current set */
    int *set;         /* [max_size] current set */
    int load, service, load_time;
    int *dp;        /* [2^max_size * max_size] Held-Karp table */
    int8_t *parent; /* [2^max_size * max_size] predecessor in the optimal path */
    int *tour;      /* [max_size] optimal order of the current set */
} Builder;

/* Optimal closed tour depot -> set -> depot (Held-Karp). Writes the order to b->tour. */
static int solve_tsp(Builder *b)
{
    const Instance *in = b->in;
    const int k = b->k, D = b->depot;
    const int *s = b->set;
    if (k <= 2) {
        for (int i = 0; i < k; i++)
            b->tour[i] = s[i];
        return k == 1 ? 2 * inst_dist(in, D, s[0])
                      : inst_dist(in, D, s[0]) + inst_dist(in, s[0], s[1]) + inst_dist(in, s[1], D);
    }
    const int full = (1 << k) - 1;
    int *dp = b->dp;
    int8_t *par = b->parent;
    for (int m = 0; m <= full; m++)
        for (int j = 0; j < k; j++)
            dp[m * k + j] = INT_MAX;
    for (int j = 0; j < k; j++) {
        dp[(1 << j) * k + j] = inst_dist(in, D, s[j]);
        par[(1 << j) * k + j] = -1;
    }
    for (int m = 1; m <= full; m++)
        for (int j = 0; j < k; j++) {
            const int cur = dp[m * k + j];
            if (cur == INT_MAX || !(m & (1 << j)))
                continue;
            for (int n = 0; n < k; n++) {
                if (m & (1 << n))
                    continue;
                const int mm = m | (1 << n), cand = cur + inst_dist(in, s[j], s[n]);
                if (cand < dp[mm * k + n]) {
                    dp[mm * k + n] = cand;
                    par[mm * k + n] = (int8_t)j;
                }
            }
        }
    int best = INT_MAX, last = 0;
    for (int j = 0; j < k; j++) {
        const int len = dp[full * k + j] + inst_dist(in, s[j], D);
        if (len < best) {
            best = len;
            last = j;
        }
    }
    for (int m = full, j = last, pos = k - 1; j >= 0; pos--) {
        b->tour[pos] = s[j];
        const int pj = par[m * k + j];
        m &= ~(1 << j);
        j = pj;
    }
    return best;
}

/* Appends the current set of the builder to the pool: sequences it with solve_tsp() and computes
 * its distance, durations and the vessel types that can carry it. */
static void record_route(Builder *b)
{
    RoutePool *pool = b->pool;
    const Instance *in = b->in;
    if (pool->count == pool->capacity) {
        pool->capacity = pool->capacity ? 2 * pool->capacity : 1024;
        pool->route = xrealloc(pool->route, sizeof(Route) * (size_t)pool->capacity);
        pool->stop = xrealloc(pool->stop, sizeof(int) * (size_t)(pool->capacity * pool->max_size));
    }
    Route *r = &pool->route[pool->count];
    r->size = (uint8_t)b->k;
    r->load = b->load;
    r->load_time = b->load_time;
    r->dist = solve_tsp(b);
    /* Same float evaluation as the original implementation. */
    r->time = (int)(((float)r->dist * in->type[0].time_per_dist + (float)b->service) + (float)b->load_time);
    r->span = (uint8_t)imin(r->time / SECONDS_PER_DAY, UINT8_MAX);
    r->types = 0;
    for (int t = 0; t < in->num_types; t++)
        if (b->load <= in->type[t].capacity) /* Q_r <= U_v, as in the paper and in the heuristic */
            r->types |= (uint16_t)(1u << t);
    memcpy(&pool->stop[pool->count * pool->max_size], b->tour, sizeof(int) * (size_t)b->k);
    pool->count++;
}

/* Depth-first enumeration of the subsets of members[from..n-1] extending the current set. */
static void enumerate(Builder *b, const int *members, int n, int from)
{
    const Instance *in = b->in;
    for (int i = from; i < n; i++) {
        const int u = members[i], v = in->inst[u].first_visit;
        const int load = in->visit_load[v];
        if (b->load + load > b->max_capacity)
            continue;
        b->set[b->k++] = u;
        b->load += load;
        b->service += in->visit_time[v];
        b->load_time += in->visit_load_time[v * in->num_berths + 0];
        record_route(b);
        if (b->k < b->pool->max_size)
            enumerate(b, members, n, i + 1);
        b->k--;
        b->load -= load;
        b->service -= in->visit_time[v];
        b->load_time -= in->visit_load_time[v * in->num_berths + 0];
    }
}

RoutePool *route_pool_build(const Instance *in, int max_size)
{
    const double t0 = process_cpu_time();
    RoutePool *pool = XCALLOC_ARRAY(RoutePool, 1);
    pool->inst = in;
    pool->max_size = max_size;
    if (in->num_types > 16)
        die("the exact models support at most 16 vessel types");
    for (int b = 1; b < in->num_berths; b++)
        if (in->berth[b].setup_time != in->berth[0].setup_time || in->berth[b].load_rate != in->berth[0].load_rate) {
            fprintf(stderr,
                    "warning: berths have different loading parameters; the exact models use those of "
                    "berth %d\n",
                    in->berth[0].id);
            break;
        }
    for (int t = 1; t < in->num_types; t++)
        if (in->type[t].time_per_dist != in->type[0].time_per_dist) {
            fprintf(stderr, "warning: vessel types have different speeds; the exact models use that of type %d\n",
                    in->type[0].id);
            break;
        }

    Builder b = {.pool = pool, .in = in, .depot = inst_berth_loc(in, 0)};
    for (int t = 0; t < in->num_types; t++)
        b.max_capacity = imax(b.max_capacity, in->type[t].capacity);
    b.set = XMALLOC_ARRAY(int, max_size);
    b.tour = XMALLOC_ARRAY(int, max_size);
    b.dp = XMALLOC_ARRAY(int, ((size_t)1 << max_size) * (size_t)max_size);
    b.parent = XMALLOC_ARRAY(int8_t, ((size_t)1 << max_size) * (size_t)max_size);
    b.load_time = in->berth[0].setup_time;

    int *members = XMALLOC_ARRAY(int, in->num_installations);
    for (int c = 0; c < in->num_clusters; c++) {
        int n = 0;
        for (int i = 0; i < in->num_installations; i++)
            if (in->inst[i].cluster == c)
                members[n++] = i;
        enumerate(&b, members, n, 0);
    }
    free(members);
    free(b.set);
    free(b.tour);
    free(b.dp);
    free(b.parent);
    pool->cpu_seconds = process_cpu_time() - t0;
    return pool;
}

void route_pool_free(RoutePool *pool)
{
    if (pool == nullptr)
        return;
    free(pool->route);
    free(pool->stop);
    free(pool);
}

bool route_visits(const RoutePool *pool, long r, int inst)
{
    const int *s = route_stops(pool, r);
    for (int i = 0; i < pool->route[r].size; i++)
        if (s[i] == inst)
            return true;
    return false;
}

void route_pool_print_stats(const RoutePool *pool)
{
    long *hist = XCALLOC_ARRAY(long, pool->max_size + 1);
    long long dist = 0;
    for (long r = 0; r < pool->count; r++) {
        hist[pool->route[r].size]++;
        dist += pool->route[r].dist;
    }
    printf("voyages: %ld (enumerated and sequenced in %.2f s)\n", pool->count, pool->cpu_seconds);
    printf("voyages per number of installations:");
    for (int k = 1; k <= pool->max_size; k++)
        printf(" %d:%ld", k, hist[k]);
    printf("\nsum of voyage distances: %lld\n", dist);
    free(hist);
}
