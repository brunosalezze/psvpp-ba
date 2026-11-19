// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file routes.h
 * @brief Enumeration of the voyage set \f$R\f$ of the exact models (paper, Section 2.1).
 *
 * A voyage is a set of installations of one cluster, of at most
 * ExactParams::max_route_size elements, whose total load fits in the largest
 * vessel (\f$Q_r \le U_v\f$; a voyage may fill its vessel exactly). Its sequence is an optimal travelling salesman tour
 * from the supply base, computed exactly by dynamic programming (Held-Karp), and its durations follow the paper:
 * \f$T_r + D_r / S_v\f$ for sailing and service, plus the berth loading time
 * \f$L^f_b + Q_r L^q_b\f$.
 *
 * The exact models treat the berths as identical and all vessel types as
 * equally fast (true for the benchmark instances); the parameters of the
 * first berth and the speed of the first type are used.
 */
#ifndef PSVPP_EXACT_ROUTES_H
#define PSVPP_EXACT_ROUTES_H

#include "common/instance.h"

#include <stdint.h>

/** @brief One voyage of the pool. */
typedef struct Route {
    int load;       /**< total load \f$Q_r\f$ (load units) */
    int dist;       /**< length of the optimal tour \f$D_r\f$ */
    int load_time;  /**< loading time at the berth, including the set-up time (s) */
    int time;       /**< loading + sailing + service time (s) */
    uint8_t size;   /**< number of installations */
    uint8_t span;   /**< whole days spanned after the departure day: time / day */
    uint16_t types; /**< bit k set if vessel type k can perform the voyage */
} Route;

/** @brief The set of all enumerated voyages. */
typedef struct RoutePool {
    const Instance *inst; /**< instance */
    int max_size;         /**< maximum installations per voyage */
    long count;           /**< number of voyages */
    long capacity;        /**< allocated entries */
    Route *route;         /**< [count] voyages */
    int *stop;            /**< [count * max_size] installations of each voyage, in tour order */
    double cpu_seconds;   /**< enumeration time */
} RoutePool;

/**
 * @brief Enumerates every feasible voyage of the instance.
 *
 * @param inst     instance.
 * @param max_size maximum number of installations per voyage.
 * @return The pool; release it with route_pool_free().
 */
[[nodiscard]] RoutePool *route_pool_build(const Instance *inst, int max_size);

/** @brief Releases a pool (accepts @c nullptr). */
void route_pool_free(RoutePool *pool);

/** @brief Installations of voyage @p r, in tour order (Route::size entries). */
static inline const int *route_stops(const RoutePool *pool, long r) { return &pool->stop[r * pool->max_size]; }

/** @brief True if vessel type @p type can perform voyage @p r. */
static inline bool route_fits(const Route *r, int type) { return (r->types >> type) & 1u; }

/** @brief True if voyage @p r visits installation @p inst. */
bool route_visits(const RoutePool *pool, long r, int inst);

/**
 * @brief Matrix \f$M^3\f$: is a vessel that starts voyage @p r on day
 *        @p start still busy with it on day @p day (cyclically)?
 */
static inline bool route_active(const Instance *in, const Route *r, int start, int day)
{
    int offset = day - start;
    if (offset < 0)
        offset += in->num_periods;
    return offset <= r->span;
}

/** @brief Prints the size histogram and the enumeration time of a pool. */
void route_pool_print_stats(const RoutePool *pool);

#endif /* PSVPP_EXACT_ROUTES_H */
