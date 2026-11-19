// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file instance.h
 * @brief Problem data of a PSVPP-BA instance (paper, Section 2) and its loader.
 *
 * An instance describes a set of offshore installations (\f$N^{OI}\f$), the
 * berths of the onshore supply base (\f$N^B\f$), the vessel types
 * (\f$V'\f$) and a cyclic planning horizon of @c num_periods days (\f$T\f$).
 * The text format is documented in docs/instance-format.md.
 *
 * ### Internal units (kept from the original implementation)
 * | quantity | unit                                                        |
 * |----------|-------------------------------------------------------------|
 * | time     | seconds (::SECONDS_PER_DAY per period)                      |
 * | load     | 1/1000 m^2 of deck area                                     |
 * | distance | 1000 x Euclidean distance between input coordinates         |
 * | cost     | 1e-6 of the reported cost (see ::COST_SCALE)                |
 *
 * ### Nodes
 * The heuristic works on a set of @e nodes:
 *  - @e visit nodes `0 .. num_visits-1`: installation @c i with frequency
 *    \f$f_i\f$ owns the \f$f_i\f$ consecutive nodes starting at
 *    Installation::first_visit (one node per weekly delivery);
 *  - @e head nodes `num_visits .. num_nodes-1`: one per (period, vessel slot),
 *    the depot node at which a voyage starts and ends (see inst_head_node()).
 *
 * ### Locations
 * Locations index the distance matrix: installations come first
 * (`0 .. num_installations-1`), then berths (`num_installations + b`).
 */
#ifndef PSVPP_COMMON_INSTANCE_H
#define PSVPP_COMMON_INSTANCE_H

#include <stdint.h>

/** @name Time constants and modelling constants of the paper */
///@{
enum {
    SECONDS_PER_DAY = 86400,      /**< length of one period (s) */
    MAX_DEPARTURES_PER_BERTH = 2, /**< departing positions per berth and day, \f$|K_b|\f$ */
    /** Berth occupation tolerance TOL = 1.5 periods, in seconds (paper, Section 2). */
    BERTH_TOLERANCE = SECONDS_PER_DAY + SECONDS_PER_DAY / 2,
};
///@}

/** @brief Divide an internal cost by this value to obtain the reported cost. */
#define COST_SCALE 1e6

/** @brief One vessel type \f$v \in V'\f$. */
typedef struct VesselType {
    int id;              /**< identifier from the instance file */
    int max_slots;       /**< vessels of this type that may be based at each berth (heuristic) */
    int max_vessels;     /**< upper bound on the fleet of this type (exact models) */
    int fixed_cost;      /**< weekly fixed cost \f$C^f_v\f$ (internal cost units) */
    float var_cost;      /**< cost per distance unit \f$C^d_v\f$ */
    int capacity;        /**< deck capacity \f$U_v\f$ (load units) */
    float time_per_dist; /**< sailing seconds per distance unit, \f$1/S_v\f$ */
} VesselType;

/** @brief One berth \f$b \in N^B\f$ of the supply base. */
typedef struct Berth {
    int id;          /**< identifier from the instance file */
    float x, y;      /**< coordinates */
    int setup_time;  /**< fixed manoeuvring and start-up time \f$L^f_b\f$ (s) */
    float load_rate; /**< loading seconds per load unit \f$L^q_b\f$ */
} Berth;

/** @brief One offshore installation \f$i \in N^{OI}\f$. */
typedef struct Installation {
    int id;                 /**< identifier from the instance file */
    int cluster;            /**< 0-based cluster index (routes never mix clusters) */
    float x, y;             /**< coordinates */
    int demand;             /**< weekly demand (load units); each visit carries demand / frequency */
    int service_time;       /**< fixed service time per visit (s) */
    float unload_rate;      /**< service seconds per load unit */
    int frequency;          /**< visits per week \f$f_i\f$ */
    int first_visit;        /**< node id of the first visit of this installation */
    int num_patterns;       /**< number of delivery patterns \f$|P_i|\f$ */
    const uint8_t *pattern; /**< `pattern[p * num_periods + t]` = 1 if pattern p visits on day t */
} Installation;

/** @brief Identifies a voyage: the vessel slot (berth, type, slot) and its departure day. */
typedef struct VoyageKey {
    int8_t period; /**< departure day */
    int8_t berth;  /**< berth the vessel is based at */
    int8_t slot;   /**< index of the vessel among the slots of its type at the berth */
    int8_t type;   /**< vessel type */
} VoyageKey;

/** @brief Complete instance data. Read-only after inst_load(). */
typedef struct Instance {
    char *name; /**< file name the instance was read from */

    int num_types;         /**< \f$|V'|\f$ */
    int num_berths;        /**< \f$|N^B|\f$ */
    int num_installations; /**< \f$|N^{OI}|\f$ */
    int num_periods;       /**< \f$|T|\f$ */
    int num_clusters;      /**< number of installation clusters */
    int horizon;           /**< planning horizon in seconds (num_periods days) */

    VesselType *type;      /**< [num_types] */
    Berth *berth;          /**< [num_berths] */
    Installation *inst;    /**< [num_installations] */
    uint8_t *pattern_data; /**< storage behind Installation::pattern */

    /* ---- locations and distances ---- */
    int num_locations; /**< num_installations + num_berths */
    int *dist;         /**< [num_locations^2] rounded distances */
    int max_dist;      /**< largest distance between two locations (Shaw normalisation, D) */
    int max_demand;    /**< largest weekly installation demand (Shaw normalisation, Q) */

    /* ---- vessel slots (heuristic fleet layout) ---- */
    int slots_per_berth;   /**< sum of VesselType::max_slots */
    int *type_slot_offset; /**< [num_types] first slot of each type inside a berth */
    int num_vessels;       /**< num_berths * slots_per_berth */

    /* ---- nodes ---- */
    int num_visits;       /**< total weekly visits, sum of frequencies */
    int num_nodes;        /**< num_visits + num_periods * num_vessels */
    int *node_loc;        /**< [num_nodes] location of each node */
    int *visit_inst;      /**< [num_visits] installation served by each visit node */
    int *visit_load;      /**< [num_visits] load delivered by the visit */
    int *visit_time;      /**< [num_visits] service time of the visit at the installation (s) */
    int *visit_load_time; /**< [num_visits * num_berths] loading time of the visit cargo at each berth */
    VoyageKey *head_key;  /**< [num_periods * num_vessels] voyage identified by each head node */

    /* ---- neighbourhoods (heuristic) ---- */
    int *neighbors;     /**< [num_installations * (num_installations-1)] others sorted by distance */
    int *num_neighbors; /**< [num_installations] neighbours considered by the insertion operators */
} Instance;

/**
 * @brief Reads an instance file and derives all auxiliary data.
 *
 * The derived data are the visit nodes, the head nodes of every vessel slot,
 * the rounded distance matrix, the normalisation constants of the Shaw
 * relatedness and the sorted neighbour lists of every installation.
 *
 * @param path path of the instance file (format: docs/instance-format.md).
 * @return Newly allocated instance; release it with inst_free(). Malformed
 *         files abort the program with a message pointing to the problem.
 */
[[nodiscard]] Instance *inst_load(const char *path);

/** @brief Releases an instance created by inst_load() (accepts @c nullptr). */
void inst_free(Instance *inst);

/** @brief Prints a human-readable summary of the instance. */
void inst_print_summary(const Instance *inst);

/** @brief Distance between two locations. */
static inline int inst_dist(const Instance *in, int a, int b) { return in->dist[a * in->num_locations + b]; }

/** @brief Distance between the locations of two nodes. */
static inline int inst_node_dist(const Instance *in, int u, int v)
{
    return in->dist[in->node_loc[u] * in->num_locations + in->node_loc[v]];
}

/** @brief Location index of berth @p b. */
static inline int inst_berth_loc(const Instance *in, int b) { return in->num_installations + b; }

/** @brief True if @p node is a voyage head (depot) node. */
static inline bool inst_is_head(const Instance *in, int node) { return node >= in->num_visits; }

/** @brief Index in [0, num_vessels) of a vessel slot. */
static inline int inst_vessel(const Instance *in, int berth, int type, int slot)
{
    return berth * in->slots_per_berth + in->type_slot_offset[type] + slot;
}

/** @brief Head node of the voyage of vessel slot @p vessel departing on @p period. */
static inline int inst_head_node(const Instance *in, int period, int vessel)
{
    return in->num_visits + period * in->num_vessels + vessel;
}

/** @brief Delivery-pattern bit: does pattern @p p of installation @p i visit on day @p t? */
static inline bool inst_pattern_visits(const Installation *i, int num_periods, int p, int t)
{
    return i->pattern[p * num_periods + t] != 0;
}

#endif /* PSVPP_COMMON_INSTANCE_H */
