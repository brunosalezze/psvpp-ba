// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file instance.c
 * @brief Instance loader and derived data (see instance.h).
 */
#include "common/instance.h"

#include "common/util.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Fraction of the largest distance below which installations count as neighbours. */
static const double NEIGHBOR_DISTANCE_FRACTION = 0.35;
/** Minimum number of neighbours of every installation (when that many exist). */
enum { MIN_NEIGHBORS = 50 };

/** @name Largest input values accepted, so that the internal integer units cannot overflow */
///@{
static const double MAX_FIXED_COST = 1000.0; /**< reported cost units (1e6 internal units each) */
static const double MAX_AREA = 1e5;          /**< m^2: demands and capacities */
static const double MAX_DAYS = 100.0;        /**< days: set-up, service and loading times and rates */
///@}

/* ------------------------------------------------------------------------- */
/*  Tokenizer with line tracking, so that parse errors point to the input    */
/* ------------------------------------------------------------------------- */

typedef struct Reader {
    const char *path;
    char *text;   /* whole file */
    char *pos;    /* next unread character */
    int line;     /* line of the next unread character */
    int tok_line; /* line of the last token returned */
} Reader;

/* Reads the whole file @p path into memory and prepares @p r to tokenize it; aborts if it cannot be read. */
static void reader_open(Reader *r, const char *path)
{
    FILE *f = fopen(path, "rb");
    if (f == nullptr)
        die("cannot open instance file '%s': %s", path, strerror(errno));
    size_t cap = 1 << 14, len = 0;
    char *buf = xmalloc(cap);
    size_t got;
    while ((got = fread(buf + len, 1, cap - len - 1, f)) > 0) {
        len += got;
        if (cap - len - 1 == 0)
            buf = xrealloc(buf, cap *= 2);
    }
    if (ferror(f))
        die("error while reading '%s'", path);
    fclose(f);
    buf[len] = '\0';
    *r = (Reader){.path = path, .text = buf, .pos = buf, .line = 1};
}

static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

/* Returns the next whitespace-separated token, NUL-terminated in place. */
static char *reader_token(Reader *r, const char *what)
{
    while (is_space(*r->pos)) {
        if (*r->pos == '\n')
            r->line++;
        r->pos++;
    }
    if (*r->pos == '\0')
        die("%s: unexpected end of file while reading %s", r->path, what);
    char *start = r->pos;
    r->tok_line = r->line;
    while (*r->pos != '\0' && !is_space(*r->pos))
        r->pos++;
    if (*r->pos != '\0') { /* consume the separator */
        if (*r->pos == '\n')
            r->line++;
        *r->pos++ = '\0';
    }
    return start;
}

/* Next token of @p r as a real number; aborts with file:line and @p what if it is not one. */
static double reader_double(Reader *r, const char *what)
{
    char *tok = reader_token(r, what);
    char *end;
    errno = 0;
    const double v = strtod(tok, &end);
    if (*end != '\0' || errno != 0)
        die("%s:%d: expected a number for %s, found '%s'", r->path, r->tok_line, what, tok);
    return v;
}

/* Next token of @p r as a real number in [0, @p max_value]; aborts with file:line and @p what otherwise.
 * The limit keeps the value, once converted to the internal integer units, inside an int. */
static double reader_amount(Reader *r, const char *what, double max_value)
{
    const double v = reader_double(r, what);
    if (!(v >= 0 && v <= max_value))
        die("%s:%d: %s must be between 0 and %g, found %g", r->path, r->tok_line, what, max_value, v);
    return v;
}

/* Next token of @p r as an integer >= @p min_value; aborts with file:line and @p what otherwise. */
static int reader_int(Reader *r, const char *what, int min_value)
{
    char *tok = reader_token(r, what);
    char *end;
    errno = 0;
    const long v = strtol(tok, &end, 10);
    if (*end != '\0' || errno != 0 || v < min_value || v > INT32_MAX)
        die("%s:%d: expected an integer >= %d for %s, found '%s'", r->path, r->tok_line, min_value, what, tok);
    return (int)v;
}

/* ------------------------------------------------------------------------- */
/*  Derived data                                                             */
/* ------------------------------------------------------------------------- */

/* Rounded distance used throughout the original code: 1000 * Euclidean distance,
 * rounded half down. The arithmetic is done in float exactly as originally. */
static int rounded_distance(float x1, float y1, float x2, float y2)
{
    const float dx = x1 - x2, dy = y1 - y2;
    const double val = 1000.0 * sqrt((double)(dx * dx + dy * dy));
    int d = (int)val;
    if (val - (double)d > 0.5)
        d++;
    return d;
}

/* Fills Instance::dist for every pair of locations (installations, then berths) and Instance::max_dist. */
static void build_distances(Instance *in)
{
    const int n = in->num_locations;
    in->dist = XCALLOC_ARRAY(int, n * n);
    float *x = XMALLOC_ARRAY(float, n), *y = XMALLOC_ARRAY(float, n);
    for (int i = 0; i < in->num_installations; i++) {
        x[i] = in->inst[i].x;
        y[i] = in->inst[i].y;
    }
    for (int b = 0; b < in->num_berths; b++) {
        x[inst_berth_loc(in, b)] = in->berth[b].x;
        y[inst_berth_loc(in, b)] = in->berth[b].y;
    }
    in->max_dist = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            const int d = rounded_distance(x[i], y[i], x[j], y[j]);
            in->dist[i * n + j] = in->dist[j * n + i] = d;
            in->max_dist = imax(in->max_dist, d);
        }
    free(x);
    free(y);
}

/* Creates one visit node per weekly delivery and one head node per (day, vessel slot). */
static void build_nodes(Instance *in)
{
    const int P = in->num_periods, B = in->num_berths, T = in->num_types;

    in->type_slot_offset = XMALLOC_ARRAY(int, T);
    in->slots_per_berth = 0;
    for (int k = 0; k < T; k++) {
        in->type_slot_offset[k] = in->slots_per_berth;
        in->slots_per_berth += in->type[k].max_slots;
    }
    if (in->slots_per_berth == 0)
        die("%s: no vessel type has slots at the berths (second field of the vessel types)", in->name);
    in->num_vessels = B * in->slots_per_berth;

    in->num_visits = 0;
    for (int i = 0; i < in->num_installations; i++)
        in->num_visits += in->inst[i].frequency;
    in->num_nodes = in->num_visits + P * in->num_vessels;

    in->node_loc = XMALLOC_ARRAY(int, in->num_nodes);
    in->visit_inst = XMALLOC_ARRAY(int, in->num_visits);
    in->visit_load = XMALLOC_ARRAY(int, in->num_visits);
    in->visit_time = XMALLOC_ARRAY(int, in->num_visits);
    in->visit_load_time = XMALLOC_ARRAY(int, in->num_visits * B);
    in->head_key = XMALLOC_ARRAY(VoyageKey, P * in->num_vessels);

    int node = 0;
    for (int i = 0; i < in->num_installations; i++) {
        Installation *u = &in->inst[i];
        u->first_visit = node;
        for (int k = 0; k < u->frequency; k++, node++) {
            in->node_loc[node] = i;
            in->visit_inst[node] = i;
            /* The expressions (and their float evaluation) follow the original code. */
            in->visit_load[node] = u->demand / u->frequency;
            in->visit_time[node] =
                (int)((float)u->service_time + (float)u->demand * u->unload_rate / (float)u->frequency);
            for (int b = 0; b < B; b++)
                in->visit_load_time[node * B + b] = (int)((float)in->visit_load[node] * in->berth[b].load_rate);
        }
    }
    for (int p = 0; p < P; p++)
        for (int b = 0; b < B; b++)
            for (int k = 0; k < T; k++)
                for (int l = 0; l < in->type[k].max_slots; l++) {
                    const int v = inst_vessel(in, b, k, l);
                    const int h = inst_head_node(in, p, v);
                    in->node_loc[h] = inst_berth_loc(in, b);
                    in->head_key[h - in->num_visits] =
                        (VoyageKey){.period = (int8_t)p, .berth = (int8_t)b, .slot = (int8_t)l, .type = (int8_t)k};
                }
}

/* qsort comparator context (qsort has no user pointer in ISO C). */
static const int *g_sort_dist_row;

/* qsort comparator: orders installation indices by their distance in the row g_sort_dist_row. */
static int compare_by_distance(const void *a, const void *b)
{
    const int ia = *(const int *)a, ib = *(const int *)b;
    const int da = g_sort_dist_row[ia], db = g_sort_dist_row[ib];
    return da != db ? (da > db) - (da < db) : (ia > ib) - (ia < ib); /* ties by index: qsort is not stable */
}

/* For every installation, sorts the other installations by distance (Instance::neighbors) and sets
 * how many of them the insertion operators consider: those closer than 35% of the largest
 * distance, but at least 50 (Instance::num_neighbors). */
static void build_neighbors(Instance *in)
{
    const int n = in->num_installations, m = n - 1;
    in->neighbors = XMALLOC_ARRAY(int, n * imax(m, 1));
    in->num_neighbors = XMALLOC_ARRAY(int, n);
    const double limit = NEIGHBOR_DISTANCE_FRACTION * in->max_dist;
    for (int i = 0; i < n; i++) {
        int *nb = &in->neighbors[i * m];
        for (int j = 0; j < m; j++)
            nb[j] = (j < i ? j : j + 1);
        g_sort_dist_row = &in->dist[i * in->num_locations];
        qsort(nb, (size_t)m, sizeof *nb, compare_by_distance);
        int count = 0;
        while (count < m && inst_dist(in, i, nb[count]) < limit)
            count++;
        in->num_neighbors[i] = imin(imax(count, MIN_NEIGHBORS), m);
    }
}

/* ------------------------------------------------------------------------- */
/*  Loader                                                                   */
/* ------------------------------------------------------------------------- */

Instance *inst_load(const char *path)
{
    Reader r;
    reader_open(&r, path);
    Instance *in = XCALLOC_ARRAY(Instance, 1);
    in->name = xstrdup(path);

    in->num_types = reader_int(&r, "the number of vessel types", 1);
    in->num_berths = reader_int(&r, "the number of berths", 1);
    in->num_installations = reader_int(&r, "the number of installations", 1);
    in->num_periods = reader_int(&r, "the number of periods", 1);
    in->num_clusters = reader_int(&r, "the number of clusters", 1);
    if (in->num_types > INT8_MAX || in->num_berths > INT8_MAX || in->num_periods > INT8_MAX)
        die("%s: at most %d vessel types, berths and periods are supported", path, INT8_MAX);
    in->horizon = in->num_periods * SECONDS_PER_DAY;

    in->type = XCALLOC_ARRAY(VesselType, in->num_types);
    for (int k = 0; k < in->num_types; k++) {
        VesselType *t = &in->type[k];
        t->id = reader_int(&r, "a vessel type id", 0);
        t->max_slots = reader_int(&r, "the vessel slots per berth", 0);
        t->max_vessels = reader_int(&r, "the maximum number of vessels", 0);
        t->fixed_cost = (int)(1e6 * reader_amount(&r, "the vessel fixed cost", MAX_FIXED_COST));
        t->var_cost = (float)(1e3 * reader_amount(&r, "the vessel variable cost", 1e3));
        t->capacity = (int)(1e3 * reader_amount(&r, "the vessel capacity", MAX_AREA));
        (void)reader_double(&r, "an unused vessel field");
        (void)reader_double(&r, "an unused vessel field");
        const double knots = reader_double(&r, "the vessel speed");
        if (!(knots > 0))
            die("%s: vessel type %d has a non-positive speed", path, t->id);
        t->time_per_dist = (float)(3.6 / (1.852 * knots));
        if (t->max_slots > INT8_MAX)
            die("%s: vessel type %d has too many slots", path, t->id);
    }

    in->berth = XCALLOC_ARRAY(Berth, in->num_berths);
    for (int b = 0; b < in->num_berths; b++) {
        Berth *s = &in->berth[b];
        s->id = reader_int(&r, "a berth id", 0);
        s->x = (float)reader_double(&r, "a berth x coordinate");
        s->y = (float)reader_double(&r, "a berth y coordinate");
        s->setup_time = (int)(SECONDS_PER_DAY * reader_amount(&r, "the berth set-up time", MAX_DAYS));
        s->load_rate = (float)(SECONDS_PER_DAY * reader_amount(&r, "the berth loading rate", MAX_DAYS) / 1e3);
    }

    in->inst = XCALLOC_ARRAY(Installation, in->num_installations);
    for (int i = 0; i < in->num_installations; i++) {
        Installation *u = &in->inst[i];
        u->id = reader_int(&r, "an installation id", 0);
        u->cluster = reader_int(&r, "an installation cluster", 1) - 1;
        u->x = (float)reader_double(&r, "an installation x coordinate");
        u->y = (float)reader_double(&r, "an installation y coordinate");
        u->demand = (int)(1e3 * reader_amount(&r, "an installation demand", MAX_AREA));
        (void)reader_double(&r, "an installation back-load (unused)");
        u->service_time = (int)(SECONDS_PER_DAY * reader_amount(&r, "an installation service time", MAX_DAYS));
        u->unload_rate = (float)(SECONDS_PER_DAY * reader_amount(&r, "an installation unloading rate", MAX_DAYS) / 1e6);
        u->frequency = reader_int(&r, "an installation frequency", 1);
        (void)reader_double(&r, "an installation time-window start (unused)");
        (void)reader_double(&r, "an installation time-window end (unused)");
        if (u->cluster >= in->num_clusters)
            die("%s: installation %d belongs to cluster %d but only %d clusters are declared", path, u->id,
                u->cluster + 1, in->num_clusters);
        if (u->frequency > in->num_periods)
            die("%s: installation %d has frequency %d > %d periods", path, u->id, u->frequency, in->num_periods);
    }

    /* Delivery patterns: groups of codes per frequency; each code is a bit mask
     * over the periods with the most significant bit corresponding to day 1. */
    const int P = in->num_periods;
    const int groups = reader_int(&r, "the number of pattern groups", 0);
    typedef struct {
        int freq, count;
        int *codes;
    } Group;
    Group *g = XCALLOC_ARRAY(Group, imax(groups, 1));
    int total_rows = 0;
    for (int k = 0; k < groups; k++) {
        g[k].freq = reader_int(&r, "a pattern group frequency", 1);
        g[k].count = reader_int(&r, "a pattern group size", 0);
        g[k].codes = XMALLOC_ARRAY(int, imax(g[k].count, 1));
        for (int l = 0; l < g[k].count; l++) {
            g[k].codes[l] = reader_int(&r, "a pattern code", 0);
            if (P < 31 && g[k].codes[l] >= (1 << P))
                die("%s: pattern code %d does not fit in %d periods", path, g[k].codes[l], P);
        }
    }
    for (int i = 0; i < in->num_installations; i++)
        for (int k = 0; k < groups; k++)
            if (g[k].freq == in->inst[i].frequency)
                total_rows += g[k].count;
    in->pattern_data = XCALLOC_ARRAY(uint8_t, imax(total_rows * P, 1));
    uint8_t *row = in->pattern_data;
    for (int i = 0; i < in->num_installations; i++) {
        Installation *u = &in->inst[i];
        u->pattern = row;
        for (int k = 0; k < groups; k++) {
            if (g[k].freq != u->frequency)
                continue;
            for (int l = 0; l < g[k].count; l++, row += P) {
                int code = g[k].codes[l], ones = 0;
                for (int t = P - 1; t >= 0; t--, code /= 2)
                    ones += row[t] = (uint8_t)(code % 2);
                if (ones != u->frequency)
                    die("%s: pattern code %d has %d visits but is listed for frequency %d", path, g[k].codes[l], ones,
                        u->frequency);
                u->num_patterns++;
            }
        }
        if (u->num_patterns == 0)
            die("%s: no delivery pattern is given for frequency %d (installation %d)", path, u->frequency, u->id);
        if (u->num_patterns > INT8_MAX)
            die("%s: more than %d delivery patterns for frequency %d", path, INT8_MAX, u->frequency);
    }
    for (int k = 0; k < groups; k++)
        free(g[k].codes);
    free(g);
    free(r.text);

    in->num_locations = in->num_installations + in->num_berths;
    in->max_demand = 0;
    for (int i = 0; i < in->num_installations; i++)
        in->max_demand = imax(in->max_demand, in->inst[i].demand);

    build_distances(in);
    build_nodes(in);
    build_neighbors(in);
    return in;
}

void inst_free(Instance *in)
{
    if (in == nullptr)
        return;
    free(in->name);
    free(in->type);
    free(in->berth);
    free(in->inst);
    free(in->pattern_data);
    free(in->dist);
    free(in->type_slot_offset);
    free(in->node_loc);
    free(in->visit_inst);
    free(in->visit_load);
    free(in->visit_time);
    free(in->visit_load_time);
    free(in->head_key);
    free(in->neighbors);
    free(in->num_neighbors);
    free(in);
}

void inst_print_summary(const Instance *in)
{
    printf("instance %s: %d installations (%d weekly visits), %d berths, %d vessel types, %d periods, %d cluster%s\n",
           in->name, in->num_installations, in->num_visits, in->num_berths, in->num_types, in->num_periods,
           in->num_clusters, in->num_clusters == 1 ? "" : "s");
    for (int k = 0; k < in->num_types; k++) {
        const VesselType *t = &in->type[k];
        printf("  vessel type %d: capacity %.0f m2, fixed cost %.2f, cost/distance %.6f, "
               "slots/berth %d, max fleet %d\n",
               t->id, t->capacity / 1e3, t->fixed_cost / COST_SCALE, t->var_cost / COST_SCALE * 1e3, t->max_slots,
               t->max_vessels);
    }
}
