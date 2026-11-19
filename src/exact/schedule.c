// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file schedule.c
 * @brief Vessel and berth scheduling of a set of voyages (see schedule.h).
 */
#include "exact/schedule.h"

#include "common/util.h"

#include <stdlib.h>
#include <string.h>

/** Nodes after which schedule_plan() gives up (it normally needs a handful). */
enum { PLAN_MAX_NODES = 5000 };
/** Days and assignments after which schedule_berths() gives up (it normally needs a few dozen). */
enum { BERTH_MAX_WORK = 200000 };
/** Longest window of consecutive days tested on its own before the whole week. */
enum { BERTH_WINDOW = 3 };

/* Returns @p buf grown to hold at least @p count elements of @p size bytes (the contents are kept). */
static void *reserve(void *buf, int *capacity, int count, size_t size)
{
    if (count <= *capacity)
        return buf;
    int cap = *capacity ? *capacity : 16;
    while (cap < count)
        cap *= 2;
    *capacity = cap;
    return xrealloc(buf, size * (size_t)cap);
}
#define RESERVE(ptr, capacity, count) ((ptr) = reserve((ptr), &(capacity), (count), sizeof *(ptr)))

/* ------------------------------------------------------------------------- */
/*  Workspace                                                                */
/* ------------------------------------------------------------------------- */

/* One berth in a state of the search over the days (times in seconds).
 *
 * The week is cyclic, so the waiting time x that the berth carries into day 0 is the one it
 * carries out of the last day. Starting the week with an idle berth instead gives the carries
 * `carry`; with the true x the carry into a day is max(carry, x + sum), where `sum` accumulates
 * loading time minus one day per day (Lindley's recursion). Every limit on the true carry of a day
 * is therefore a limit on `carry` and an upper bound `room` on x. At the end of the week x is
 * known: it is the last `carry`, and it must not exceed `room`. Once x + sum <= carry holds for
 * every admissible x, the start of the week no longer matters and `sum` is dropped (RESOLVED). */
typedef struct Chain {
    int32_t carry; /* waiting time carried into the current day when the week starts idle */
    int32_t sum;   /* loading time minus elapsed days since day 0, or RESOLVED */
    int32_t room;  /* largest carry into day 0 allowed by the days so far */
} Chain;

static const int32_t RESOLVED = INT32_MIN;

/* The states at the start of one day from which no schedule exists. */
typedef struct Failed {
    Chain *chain; /* [count * B] */
    int count, cap_chain;
} Failed;

struct Scheduler {
    ScheduleStats stats;

    /* ---- berths ---- */
    const ScheduleProblem *p;
    const int *max_wait;
    int B, K;
    int *order, *first;  /* voyages sorted by day of the rotated week; first voyage of each day */
    Chain *path;         /* [(P + 1) * B] berths at the start of each day on the current path, sorted */
    uint8_t *path_from;  /* [P * B] berth of the day that each berth of the next day continues */
    int *path_voyage;    /* [P * B * K] voyages that it loads on the day, in loading order (-1: none) */
    int *group, *sorted; /* [P * B * K] voyages given to each berth of the day so far; in loading order */
    int *size;           /* [P * B] their number */
    int *choice;         /* [n] berth chosen for each voyage */
    int64_t *key;        /* [n * B] preference of each berth for each voyage (-1: not allowed or tried) */
    Chain *next;         /* [B] berths after the day */
    int *perm, *label, *label2;
    int64_t *max_total;            /* [span + 1] largest total carry into each day that the following days can absorb */
    int64_t *day_load;             /* [P] loading time of each day of the rotated week */
    int shift;                     /* day of the week that is day 0 of the rotated week */
    int start, span;               /* days of the current search: span days from day start of the rotated week */
    bool closing;                  /* the search covers the whole week, which must close */
    int deepest;                   /* deepest day that the search had to look at */
    bool end_used;                 /* the search failed somewhere because of the days after or of the closing */
    int window_first, window_days; /* days that suffice to prove the last infeasibility */
    Failed *failed;                /* [P] */
    int num_failed_days;
    long work;      /* days and assignments examined by the current call */
    bool limit_hit; /* the current call stopped at the work limit */
    int cap_order, cap_first, cap_path, cap_path_from, cap_path_voyage, cap_group, cap_sorted, cap_size, cap_choice,
        cap_key, cap_next, cap_perm, cap_label, cap_label2, cap_max_total, cap_day_load;

    /* ---- vessels ---- */
    int periods, fleet, arcs;
    int *arc, *arc_first, *arc_end; /* arcs to place, in order: index, first day and end in the rotated week */
    int *busy, *limit, *chosen;     /* [fleet] busy until; taken from; [arcs] vessel of each arc */
    int *cover, *cross;             /* [periods] */
    uint16_t *code;                 /* [fleet + 1] key of a vessel state */
    uint16_t *dead;                 /* keys of the states without solution */
    int num_dead, key_len;
    int *bucket;
    int num_buckets;
    int cap_arc, cap_arc_first, cap_arc_end, cap_busy, cap_limit, cap_chosen, cap_cover, cap_cross, cap_code, cap_dead,
        cap_bucket;

    /* ---- plan ---- */
    int *base, *slack, *reach, *wait, *trial; /* [n] per voyage, see schedule_plan() */
    uint8_t *extra, *most;                    /* [n] extra days granted / useful */
    bool *inside;                             /* [n] voyages of the days in conflict */
    int *type_start, *type_length, *type_index, *type_vessel;
    uint8_t *seen; /* extra-day vectors already visited */
    int num_seen;
    int *seen_bucket;
    int num_seen_buckets;
    uint8_t *conflict; /* berth conflicts found: per voyage the extra days it had, or NOT_IN_CONFLICT */
    int num_conflicts;
    long nodes;
    bool undecided;
    int cap_base, cap_slack, cap_reach, cap_wait, cap_trial, cap_extra, cap_most, cap_type_start, cap_type_length,
        cap_type_index, cap_type_vessel, cap_seen, cap_seen_bucket, cap_conflict, cap_inside;
};

Scheduler *scheduler_new(void) { return XCALLOC_ARRAY(Scheduler, 1); }

void scheduler_free(Scheduler *s)
{
    if (s == nullptr)
        return;
    for (int i = 0; i < s->num_failed_days; i++)
        free(s->failed[i].chain);
    void *const buffers[] = {
        s->failed,     s->order,       s->first,      s->path,        s->path_from, s->path_voyage, s->group,
        s->sorted,     s->size,        s->choice,     s->key,         s->next,      s->perm,        s->label,
        s->label2,     s->max_total,   s->day_load,   s->arc,         s->arc_first, s->arc_end,     s->busy,
        s->limit,      s->chosen,      s->cover,      s->cross,       s->code,      s->dead,        s->bucket,
        s->base,       s->slack,       s->reach,      s->wait,        s->trial,     s->extra,       s->most,
        s->type_start, s->type_length, s->type_index, s->type_vessel, s->seen,      s->seen_bucket, s->conflict,
        s->inside};
    for (size_t i = 0; i < sizeof buffers / sizeof buffers[0]; i++)
        free(buffers[i]);
    free(s);
}

ScheduleStats scheduler_stats(const Scheduler *s) { return s->stats; }

bool scheduler_gave_up(const Scheduler *s) { return s->limit_hit; }

/* FNV-1a hash of @p len bytes. */
static uint32_t hash_bytes(const void *data, size_t len)
{
    const unsigned char *b = data;
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < len; i++)
        h = (h ^ b[i]) * 16777619u;
    return h;
}

/* ------------------------------------------------------------------------- */
/*  Berths                                                                   */
/* ------------------------------------------------------------------------- */

/* Waiting-time limit of voyage @p v. */
static int wait_limit(const Scheduler *s, int v) { return s->max_wait ? s->max_wait[v] : WAIT_UNLIMITED; }

/* Voyages that are interchangeable for the berths. */
static bool same_voyage(const Scheduler *s, int a, int b)
{
    return s->p->load[a] == s->p->load[b] && wait_limit(s, a) == wait_limit(s, b);
}

static bool same_chain(const Chain *a, const Chain *b)
{
    return a->carry == b->carry && a->sum == b->sum && a->room == b->room;
}

/* Lexicographic order of the berths of a state. */
static bool chain_less(const Chain *a, const Chain *b)
{
    if (a->carry != b->carry)
        return a->carry < b->carry;
    if (a->sum != b->sum)
        return a->sum < b->sum;
    return a->room < b->room;
}

/* A group of voyages loaded by one berth on one day: puts them in the best loading order in
 * @p sorted (by limit plus loading time, the earliest-deadline order), and returns their total
 * loading time and the largest carry into the day that keeps every wait within its limit (it may
 * be negative: no order works). */
static void group_limits(const Scheduler *s, const int *g, int size, int *sorted, int *load, int *cap)
{
    const ScheduleProblem *p = s->p;
    for (int i = 0; i < size; i++) { /* insertion sort by limit + load, unlimited voyages last */
        const int v = g[i];
        const int64_t key = wait_limit(s, v) == WAIT_UNLIMITED ? INT64_MAX : (int64_t)wait_limit(s, v) + p->load[v];
        int j = i;
        for (; j > 0; j--) {
            const int u = sorted[j - 1];
            const int64_t other =
                wait_limit(s, u) == WAIT_UNLIMITED ? INT64_MAX : (int64_t)wait_limit(s, u) + p->load[u];
            if (other <= key)
                break;
            sorted[j] = u;
        }
        sorted[j] = v;
    }
    int64_t ahead = 0, most = p->max_carry;
    for (int i = 0; i < size; i++) {
        const int v = sorted[i];
        if (wait_limit(s, v) != WAIT_UNLIMITED && wait_limit(s, v) - ahead < most)
            most = wait_limit(s, v) - ahead;
        ahead += p->load[v];
    }
    *load = (int)ahead;
    *cap = most < INT32_MIN / 2 ? INT32_MIN / 2 : (int)most;
}

/* State of a berth after a day on which it loads @p load seconds and may start with a carry of
 * at most @p cap. Returns false if the day is infeasible for the berth. */
static bool chain_step(const ScheduleProblem *p, Chain c, int load, int cap, Chain *out)
{
    if (c.carry > cap)
        return false;
    if (c.sum != RESOLVED) { /* max(carry, x + sum) <= cap */
        if ((int64_t)cap - c.sum < c.room)
            c.room = (int32_t)(cap - c.sum);
        if (c.room < 0)
            return false;
    }
    const int64_t over = (int64_t)c.carry + load - p->day_length;
    c.carry = over > 0 ? (int32_t)over : 0;
    if (c.carry > p->max_carry)
        return false;
    if (c.sum != RESOLVED) {
        c.sum += load - p->day_length;
        if ((int64_t)p->max_carry - c.sum < c.room) /* the carry into the next day is limited too */
            c.room = (int32_t)(p->max_carry - c.sum);
        if (c.room < 0)
            return false;
        if ((int64_t)c.room + c.sum <= c.carry)
            c.sum = RESOLVED;
    }
    *out = c;
    return true;
}

/* true if state @p x is at least as good as state @p y on every berth (both sorted). */
static bool dominates(const Chain *x, const Chain *y, int berths)
{
    for (int j = 0; j < berths; j++)
        if (x[j].carry > y[j].carry || x[j].sum > y[j].sum || x[j].room < y[j].room)
            return false;
    return true;
}

static bool place_voyage(Scheduler *s, int i, int t);

/* Day of the rotated week examined at depth @p i of the current search. */
static int day_at(const Scheduler *s, int i) { return (s->start + i) % s->p->periods; }

/* Searches a schedule of the day at depth @p i and of the following ones, from the state at
 * s->path[i]. States from which the search failed are remembered: a state that is no better than
 * one of them on any berth fails too. */
static bool solve_day(Scheduler *s, int i)
{
    const int B = s->B;
    if (++s->work > BERTH_MAX_WORK) {
        s->limit_hit = true;
        return false;
    }
    const Chain *state = &s->path[i * B];
    Failed *f = &s->failed[i];
    s->work += f->count / 8; /* the comparison with the failed states is part of the work */
    for (int k = 0; k < f->count; k++)
        if (dominates(&f->chain[k * B], state, B))
            return false;
    if (i > s->deepest)
        s->deepest = i;
    for (int c = 0; c < B; c++)
        s->size[i * B + c] = 0;
    if (place_voyage(s, i, 0))
        return true;
    if (s->limit_hit)
        return false;
    for (int k = 0; k < f->count;) /* the failed states that are no better than this one are superseded */
        if (dominates(state, &f->chain[k * B], B)) {
            f->count--;
            memmove(&f->chain[k * B], &f->chain[f->count * B], sizeof(Chain) * (size_t)B);
        } else {
            k++;
        }
    RESERVE(f->chain, f->cap_chain, (f->count + 1) * B);
    memcpy(&f->chain[f->count * B], state, sizeof(Chain) * (size_t)B);
    f->count++;
    return false;
}

/* Every voyage of the day at depth @p i has a berth: computes the state of the next day and
 * continues from it. On the last day of a whole week the berths must close the week. */
static bool finish_day(Scheduler *s, int i)
{
    const ScheduleProblem *p = s->p;
    const int B = s->B, K = s->K;
    const Chain *chain = &s->path[i * B];
    const int *size = &s->size[i * B];
    int *sorted = &s->sorted[i * B * K];
    int64_t total = 0;
    s->work++;
    for (int c = 0; c < B; c++) {
        int load, cap;
        group_limits(s, &s->group[(i * B + c) * K], size[c], &sorted[c * K], &load, &cap);
        if (!chain_step(p, chain[c], load, cap, &s->next[c]))
            return false;
        total += s->next[c].carry;
    }
    if (total > s->max_total[i + 1]) {
        s->end_used = true;
        return false; /* more waiting time than the berths can absorb in the days that follow */
    }
    s->stats.berth_states++;
    const bool last = i == s->span - 1;
    if (last) {
        for (int c = 0; c < B; c++) {
            if (s->closing &&
                (s->next[c].carry > s->next[c].room || (s->next[c].sum != RESOLVED && s->next[c].sum > 0))) {
                s->end_used = true;
                return false; /* the berth cannot start the week with the carry it ends with */
            }
            s->perm[c] = c;
        }
    } else {
        for (int c = 0; c < B; c++) { /* insertion sort of the berths */
            int j = c;
            for (; j > 0 && chain_less(&s->next[c], &s->next[s->perm[j - 1]]); j--)
                s->perm[j] = s->perm[j - 1];
            s->perm[j] = c;
        }
    }
    for (int j = 0; j < B; j++) {
        const int c = s->perm[j];
        s->path[(i + 1) * B + j] = s->next[c];
        s->path_from[i * B + j] = (uint8_t)c;
        for (int k = 0; k < K; k++)
            s->path_voyage[(i * B + j) * K + k] = k < size[c] ? sorted[c * K + k] : -1;
    }
    return last || solve_day(s, i + 1);
}

/* Gives a berth to voyage number @p t of the day at depth @p i and to the following ones, in every
 * way that differs by more than an exchange of identical berths or of identical voyages. The
 * berths are tried by increasing carry left to the next day, then by increasing occupation. */
static bool place_voyage(Scheduler *s, int i, int t)
{
    const int B = s->B, K = s->K, d = day_at(s, i);
    if (t == s->first[d + 1] - s->first[d])
        return finish_day(s, i);
    const int at = s->first[d] + t, v = s->order[at];
    const Chain *chain = &s->path[i * B];
    int *group = &s->group[i * B * K], *size = &s->size[i * B], *sorted = &s->sorted[i * B * K];
    int64_t *key = &s->key[(size_t)at * (size_t)B];
    const bool twin = t > 0 && same_voyage(s, s->order[at - 1], v);
    for (int c = 0; c < B; c++) {
        key[c] = -1;
        if (size[c] == K || (twin && c < s->choice[at - 1]))
            continue;
        if (size[c] == 0 && c > 0 && size[c - 1] == 0 && same_chain(&chain[c - 1], &chain[c]))
            continue; /* an identical empty berth comes first */
        group[c * K + size[c]] = v;
        int load, cap;
        Chain after;
        group_limits(s, &group[c * K], size[c] + 1, &sorted[c * K], &load, &cap);
        if (chain_step(s->p, chain[c], load, cap, &after))
            key[c] = ((int64_t)after.carry << 32) + chain[c].carry + load - s->p->load[v];
    }
    for (;;) {
        int best = -1;
        for (int c = 0; c < B; c++)
            if (key[c] >= 0 && (best < 0 || key[c] < key[best]))
                best = c;
        if (best < 0)
            return false;
        key[best] = -1;
        group[best * K + size[best]++] = v;
        s->choice[at] = best;
        if (place_voyage(s, i, t + 1))
            return true;
        size[best]--;
        if (s->limit_hit)
            return false;
    }
}

/* Searches a schedule of @p span consecutive days from day @p start of the rotated week. A whole
 * week (@p closing) must close: @p unknown berths may carry something from its last day into its
 * first, the others not. A shorter window starts with idle berths and ends freely. */
static bool search_days(Scheduler *s, int start, int span, bool closing, int unknown)
{
    const ScheduleProblem *p = s->p;
    const int B = s->B;
    s->start = start;
    s->span = span;
    s->closing = closing;
    s->deepest = 0;
    s->end_used = false;
    /* total carry that the following days can absorb (each berth unloads at most a day per day) */
    s->max_total[span] = closing ? (int64_t)unknown * p->max_carry : (int64_t)B * p->max_carry;
    for (int i = span - 1; i >= 0; --i) {
        const int64_t room = s->max_total[i + 1] + (int64_t)B * p->day_length - s->day_load[day_at(s, i)];
        s->max_total[i] = room < (int64_t)B * p->max_carry ? room : (int64_t)B * p->max_carry;
    }
    if (s->max_total[0] < 0) {
        s->end_used = true;
        s->deepest = span - 1;
        return false;
    }
    for (int i = 0; i < span; i++)
        s->failed[i].count = 0;
    for (int j = 0; j < B; j++)
        s->path[j] = j < B - unknown ? (Chain){.carry = 0, .sum = RESOLVED, .room = closing ? 0 : p->max_carry}
                                     : (Chain){.carry = 0, .sum = 0, .room = p->max_carry};
    return solve_day(s, 0);
}

/* Records the days whose voyages suffice to prove that no schedule exists. */
static void set_window(Scheduler *s, int start, int span)
{
    s->window_first = (start + s->shift) % s->p->periods;
    s->window_days = span;
}

bool schedule_berths(Scheduler *s, const ScheduleProblem *p, const bool *include, const int *max_wait, int *berth,
                     int *position)
{
    const int P = p->periods, B = p->berths, K = p->positions, n = p->n;
    s->stats.berth_calls++;
    if (B > UINT8_MAX)
        die("schedule_berths: more than %d berths are not supported", UINT8_MAX);
    s->p = p;
    s->max_wait = max_wait;
    s->B = B;
    s->K = K;
    s->work = 0;
    s->limit_hit = false;
    s->shift = 0;
    set_window(s, 0, P);

    /* voyages per day; the lightest day becomes the last day of the rotated week */
    RESERVE(s->first, s->cap_first, P + 1);
    RESERVE(s->cover, s->cap_cover, P);
    RESERVE(s->cross, s->cap_cross, P + 1);
    RESERVE(s->order, s->cap_order, imax(n, 1));
    RESERVE(s->day_load, s->cap_day_load, P);
    int *count = s->cross, *total = s->cover;
    for (int d = 0; d < P; d++)
        count[d] = total[d] = 0;
    int scheduled = 0;
    for (int v = 0; v < n; v++)
        if (include == nullptr || include[v]) {
            count[p->day[v]]++;
            total[p->day[v]] += p->load[v];
            scheduled++;
        }
    if (scheduled == 0)
        return true;
    int lightest = 0;
    for (int d = 0; d < P; d++) {
        if (count[d] > K * B) {
            set_window(s, d, 1);
            return false;
        }
        if (count[d] < count[lightest] || (count[d] == count[lightest] && total[d] < total[lightest]))
            lightest = d;
    }
    /* A berth without voyage on the last day ends the week idle (the carry is less than a day), so at
     * most `open` berths carry something into day 0. */
    const int open = p->max_carry < p->day_length ? imin(B, count[lightest]) : B;
    s->shift = (lightest + 1) % P; /* day `shift` becomes day 0 of the rotated week */
    for (int d = 0; d < P; d++)
        s->day_load[d] = total[(d + s->shift) % P];
    s->first[0] = 0; /* first[d]: voyages on the rotated days before d */
    for (int d = 0; d < P; d++)
        s->first[d + 1] = s->first[d] + count[(d + s->shift) % P];
    for (int d = 0; d < P; d++)
        count[d] = s->first[d]; /* from here on: next free index of each rotated day */
    for (int v = 0; v < n; v++)
        if (include == nullptr || include[v])
            s->order[count[(p->day[v] - s->shift + P) % P]++] = v;
    for (int d = 0; d < P; d++) /* each day sorted: heaviest first, identical voyages together */
        for (int i = s->first[d] + 1; i < s->first[d + 1]; i++) {
            const int v = s->order[i];
            int j = i;
            for (; j > s->first[d]; j--) {
                const int u = s->order[j - 1];
                if (p->load[u] > p->load[v] || (p->load[u] == p->load[v] && wait_limit(s, u) <= wait_limit(s, v)))
                    break;
                s->order[j] = u;
            }
            s->order[j] = v;
        }

    /* working memory */
    if (s->num_failed_days < P) {
        s->failed = xrealloc(s->failed, sizeof(Failed) * (size_t)P);
        memset(&s->failed[s->num_failed_days], 0, sizeof(Failed) * (size_t)(P - s->num_failed_days));
        s->num_failed_days = P;
    }
    RESERVE(s->path, s->cap_path, (P + 1) * B);
    RESERVE(s->path_from, s->cap_path_from, P * B);
    RESERVE(s->path_voyage, s->cap_path_voyage, P * B * K);
    RESERVE(s->group, s->cap_group, P * B * K);
    RESERVE(s->sorted, s->cap_sorted, P * B * K);
    RESERVE(s->size, s->cap_size, P * B);
    RESERVE(s->choice, s->cap_choice, scheduled);
    RESERVE(s->key, s->cap_key, scheduled * B);
    RESERVE(s->next, s->cap_next, B);
    RESERVE(s->perm, s->cap_perm, B);
    RESERVE(s->label, s->cap_label, B);
    RESERVE(s->label2, s->cap_label2, B);
    RESERVE(s->max_total, s->cap_max_total, P + 1);

    /* Short windows first: a few consecutive days must be loadable by berths that start idle. This
     * settles most infeasible cases quickly and tells which days are in conflict. A window that
     * begins or ends with an empty day is covered by a shorter one. */
    for (int span = 1; span <= BERTH_WINDOW && span < P; span++)
        for (int start = 0; start < P; start++) {
            const int end = (start + span - 1) % P;
            if (s->first[start + 1] == s->first[start] || s->first[end + 1] == s->first[end])
                continue;
            if (search_days(s, start, span, false, 0))
                continue;
            if (s->limit_hit)
                s->stats.berth_limits++;
            else
                set_window(s, start, span);
            return false;
        }

    /* The whole week. The carry of a berth into day 0 is the one it ends the week with: unknown at
     * this point for `open` berths, zero for the others. A first pass assumes zero for every berth,
     * which is the usual case and has far fewer states. The second pass is the complete one; it is
     * needed only if the first failed because of the end of the week. */
    bool found = search_days(s, 0, P, true, 0);
    if (!found && !s->limit_hit) {
        if (!s->end_used)
            set_window(s, 0, s->deepest + 1);
        else if (open > 0)
            found = search_days(s, 0, P, true, open);
    }
    if (s->limit_hit)
        s->stats.berth_limits++;
    if (!found)
        return false;

    if (berth || position) { /* follow the berths back through the week */
        for (int j = 0; j < B; j++)
            s->label[j] = j;
        for (int d = P - 1; d >= 0; --d) {
            for (int j = 0; j < B; j++) {
                for (int i = 0; i < K; i++) {
                    const int v = s->path_voyage[(d * B + j) * K + i];
                    if (v < 0)
                        continue;
                    if (berth)
                        berth[v] = s->label[j];
                    if (position)
                        position[v] = i;
                }
                s->label2[s->path_from[d * B + j]] = s->label[j];
            }
            int *swap = s->label;
            s->label = s->label2;
            s->label2 = swap;
            const int cap = s->cap_label;
            s->cap_label = s->cap_label2;
            s->cap_label2 = cap;
        }
    }
    return true;
}

void scheduler_conflict_days(const Scheduler *s, int *first, int *days)
{
    *first = s->window_first;
    *days = s->window_days;
}

/* ------------------------------------------------------------------------- */
/*  Vessels                                                                  */
/* ------------------------------------------------------------------------- */

/* Looks up the key in s->code among the failed states; adds it when @p add is set. */
static bool failed_state(Scheduler *s, bool add)
{
    const int len = s->key_len;
    const size_t bytes = sizeof(uint16_t) * (size_t)len;
    if (s->num_buckets == 0) {
        s->num_buckets = 256;
        RESERVE(s->bucket, s->cap_bucket, s->num_buckets);
        for (int i = 0; i < s->num_buckets; i++)
            s->bucket[i] = -1;
    }
    uint32_t h = hash_bytes(s->code, bytes) & (uint32_t)(s->num_buckets - 1);
    for (; s->bucket[h] >= 0; h = (h + 1) & (uint32_t)(s->num_buckets - 1))
        if (memcmp(&s->dead[(size_t)s->bucket[h] * (size_t)len], s->code, bytes) == 0)
            return true;
    if (!add)
        return false;
    RESERVE(s->dead, s->cap_dead, (s->num_dead + 1) * len);
    memcpy(&s->dead[(size_t)s->num_dead * (size_t)len], s->code, bytes);
    s->bucket[h] = s->num_dead++;
    if (2 * s->num_dead > s->num_buckets) { /* rehash into a table twice as large */
        s->num_buckets *= 2;
        RESERVE(s->bucket, s->cap_bucket, s->num_buckets);
        for (int i = 0; i < s->num_buckets; i++)
            s->bucket[i] = -1;
        for (int k = 0; k < s->num_dead; k++) {
            uint32_t g = hash_bytes(&s->dead[(size_t)k * (size_t)len], bytes) & (uint32_t)(s->num_buckets - 1);
            while (s->bucket[g] >= 0)
                g = (g + 1) & (uint32_t)(s->num_buckets - 1);
            s->bucket[g] = k;
        }
    }
    return false;
}

/* Key of the vessel states before arc @p i: the arc number and, sorted, what every vessel can
 * still do (the vessels themselves are interchangeable). */
static void vessel_key(Scheduler *s, int i)
{
    const int from = s->arc_first[i];
    s->code[0] = (uint16_t)i;
    for (int v = 0; v < s->fleet; v++) {
        const uint16_t code = (uint16_t)(imax(s->busy[v], from) * (s->periods + 1) + s->limit[v]);
        int j = v + 1;
        for (; j > 1 && s->code[j - 1] > code; j--)
            s->code[j] = s->code[j - 1];
        s->code[j] = code;
    }
}

/* Assigns arc @p i and the following ones. */
static bool place_arc(Scheduler *s, int i)
{
    if (i == s->arcs)
        return true;
    s->stats.vessel_nodes++;
    vessel_key(s, i);
    if (failed_state(s, false))
        return false;
    const int from = s->arc_first[i], end = s->arc_end[i];
    /* free vessels that are not taken before the arc ends; one per distinct limit, tightest first */
    for (int last = -1;;) {
        int best = -1;
        for (int v = 0; v < s->fleet; v++)
            if (s->busy[v] <= from && s->limit[v] >= end && s->limit[v] > last &&
                (best < 0 || s->limit[v] < s->limit[best]))
                best = v;
        if (best < 0)
            break;
        last = s->limit[best];
        const int before = s->busy[best];
        s->busy[best] = end;
        s->chosen[i] = best;
        if (place_arc(s, i + 1))
            return true;
        s->busy[best] = before;
    }
    vessel_key(s, i);
    failed_state(s, true);
    return false;
}

bool schedule_vessels(Scheduler *s, int periods, int n, const int *start, const int *length, int fleet, int *vessel)
{
    const int P = periods;
    s->stats.vessel_calls++;
    if (n == 0)
        return true;
    if (fleet <= 0)
        return false;
    if ((int64_t)(P + 1) * (P + 1) > UINT16_MAX || n > UINT16_MAX)
        die("schedule_vessels: horizons of more than %d days are not supported", 254);
    RESERVE(s->cover, s->cap_cover, P);
    RESERVE(s->cross, s->cap_cross, P + 1);
    for (int d = 0; d < P; d++)
        s->cover[d] = s->cross[d] = 0;
    for (int a = 0; a < n; a++) {
        if (length[a] < 1 || length[a] > P)
            return false;
        for (int k = 0; k < length[a]; k++) {
            const int d = (start[a] + k) % P;
            s->cover[d]++;
            if (k > 0)
                s->cross[d]++; /* the arc contains the boundary before day d */
        }
    }
    int cut = 0;
    for (int d = 0; d < P; d++) {
        if (s->cover[d] > fleet)
            return false;
        if (s->cross[d] < s->cross[cut])
            cut = d;
    }

    /* Rotated week starting at the cut. The arcs across the cut overlap each other: they take the
     * first vessels, which are busy until the arc ends and taken again from its first day. */
    RESERVE(s->busy, s->cap_busy, fleet);
    RESERVE(s->limit, s->cap_limit, fleet);
    RESERVE(s->arc, s->cap_arc, n);
    RESERVE(s->arc_first, s->cap_arc_first, n);
    RESERVE(s->arc_end, s->cap_arc_end, n);
    RESERVE(s->chosen, s->cap_chosen, n);
    RESERVE(s->code, s->cap_code, fleet + 1);
    for (int v = 0; v < fleet; v++) {
        s->busy[v] = 0;
        s->limit[v] = P;
    }
    int wrapped = 0;
    s->arcs = 0;
    for (int a = 0; a < n; a++) {
        const int first = (start[a] - cut + P) % P;
        if (first + length[a] > P) {
            if (wrapped == fleet)
                return false;
            s->busy[wrapped] = first + length[a] - P;
            s->limit[wrapped] = first;
            if (vessel)
                vessel[a] = wrapped;
            wrapped++;
            continue;
        }
        /* insertion sort by first day, longest first */
        int j = s->arcs++;
        for (; j > 0 &&
               (s->arc_first[j - 1] > first || (s->arc_first[j - 1] == first && s->arc_end[j - 1] < first + length[a]));
             j--) {
            s->arc[j] = s->arc[j - 1];
            s->arc_first[j] = s->arc_first[j - 1];
            s->arc_end[j] = s->arc_end[j - 1];
        }
        s->arc[j] = a;
        s->arc_first[j] = first;
        s->arc_end[j] = first + length[a];
    }
    s->periods = P;
    s->fleet = fleet;
    s->key_len = fleet + 1;
    s->num_dead = 0;
    s->num_buckets = 0;
    if (!place_arc(s, 0))
        return false;
    if (vessel)
        for (int i = 0; i < s->arcs; i++)
            vessel[s->arc[i]] = s->chosen[i];
    return true;
}

/* ------------------------------------------------------------------------- */
/*  Vessels and berths together                                              */
/* ------------------------------------------------------------------------- */

enum { NOT_IN_CONFLICT = UINT8_MAX };

/* Can the vessels of @p type perform their voyages when voyage v blocks its vessel for
 * base[v] + extra[v] days? Stores the vessels in @p vessel when it is not null. */
static bool vessels_fit(Scheduler *s, const ScheduleProblem *p, int type, int *vessel)
{
    int m = 0;
    for (int v = 0; v < p->n; v++)
        if (p->type[v] == type) {
            s->type_start[m] = p->day[v];
            s->type_length[m] = s->base[v] + s->extra[v];
            s->type_index[m++] = v;
        }
    if (!schedule_vessels(s, p->periods, m, s->type_start, s->type_length, p->fleet[type],
                          vessel ? s->type_vessel : nullptr))
        return false;
    if (vessel)
        for (int i = 0; i < m; i++)
            vessel[s->type_index[i]] = s->type_vessel[i];
    return true;
}

/* Looks up the current extra-day vector among the visited ones and adds it. */
static bool seen_before(Scheduler *s, int n)
{
    if (s->num_seen_buckets == 0) {
        s->num_seen_buckets = 256;
        RESERVE(s->seen_bucket, s->cap_seen_bucket, s->num_seen_buckets);
        for (int i = 0; i < s->num_seen_buckets; i++)
            s->seen_bucket[i] = -1;
    }
    uint32_t h = hash_bytes(s->extra, (size_t)n) & (uint32_t)(s->num_seen_buckets - 1);
    for (; s->seen_bucket[h] >= 0; h = (h + 1) & (uint32_t)(s->num_seen_buckets - 1))
        if (memcmp(&s->seen[(size_t)s->seen_bucket[h] * (size_t)n], s->extra, (size_t)n) == 0)
            return true;
    RESERVE(s->seen, s->cap_seen, (s->num_seen + 1) * n);
    memcpy(&s->seen[(size_t)s->num_seen * (size_t)n], s->extra, (size_t)n);
    s->seen_bucket[h] = s->num_seen++;
    if (2 * s->num_seen > s->num_seen_buckets) {
        s->num_seen_buckets *= 2;
        RESERVE(s->seen_bucket, s->cap_seen_bucket, s->num_seen_buckets);
        for (int i = 0; i < s->num_seen_buckets; i++)
            s->seen_bucket[i] = -1;
        for (int k = 0; k < s->num_seen; k++) {
            uint32_t g = hash_bytes(&s->seen[(size_t)k * (size_t)n], (size_t)n) & (uint32_t)(s->num_seen_buckets - 1);
            while (s->seen_bucket[g] >= 0)
                g = (g + 1) & (uint32_t)(s->num_seen_buckets - 1);
            s->seen_bucket[g] = k;
        }
    }
    return false;
}

/* Waiting-time limit of voyage @p v with its current extra days. */
static int limit_with_extra(const Scheduler *s, const ScheduleProblem *p, int v)
{
    const int64_t limit = (int64_t)s->slack[v] + (int64_t)s->extra[v] * p->day_length;
    return limit >= s->reach[v] ? WAIT_UNLIMITED : (int)limit;
}

/* Searches from the current extra days (the vessels are known to fit). On success the berths are in
 * @p berth and @p position. */
static bool plan_search(Scheduler *s, const ScheduleProblem *p, int *berth, int *position)
{
    const int n = p->n;
    if (++s->nodes > PLAN_MAX_NODES) {
        s->undecided = true;
        return false;
    }
    s->stats.plan_nodes++;
    /* a berth conflict found earlier that still holds: none of its voyages has more days now */
    int hit = -1;
    for (int k = s->num_conflicts - 1; k >= 0 && hit < 0; k--) {
        const uint8_t *c = &s->conflict[(size_t)k * (size_t)n];
        bool holds = true;
        for (int v = 0; v < n && holds; v++)
            holds = c[v] == NOT_IN_CONFLICT || s->extra[v] <= c[v];
        if (holds)
            hit = k;
    }
    if (hit < 0) {
        for (int v = 0; v < n; v++)
            s->wait[v] = limit_with_extra(s, p, v);
        if (schedule_berths(s, p, nullptr, s->wait, berth, position))
            return true;
        if (s->limit_hit) {
            s->undecided = true;
            return false;
        }
        /* The voyages of a few consecutive days fail on their own. Among them, a minimal set of
         * voyages whose limits cause the failure: lift the limits one at a time and keep those
         * without which a schedule of these days appears. */
        int first, days;
        scheduler_conflict_days(s, &first, &days);
        for (int v = 0; v < n; v++)
            s->inside[v] = (p->day[v] - first + p->periods) % p->periods < days;
        memcpy(s->trial, s->wait, sizeof(int) * (size_t)n);
        for (int v = 0; v < n; v++) {
            if (!s->inside[v]) {
                s->trial[v] = WAIT_UNLIMITED;
                continue;
            }
            if (s->trial[v] == WAIT_UNLIMITED)
                continue;
            const int limit = s->trial[v];
            s->trial[v] = WAIT_UNLIMITED;
            if (schedule_berths(s, p, s->inside, s->trial, nullptr, nullptr) || s->limit_hit)
                s->trial[v] = limit; /* needed (or not proven superfluous) */
        }
        hit = s->num_conflicts;
        RESERVE(s->conflict, s->cap_conflict, (hit + 1) * n);
        for (int v = 0; v < n; v++)
            s->conflict[(size_t)hit * (size_t)n + (size_t)v] =
                s->trial[v] == WAIT_UNLIMITED ? NOT_IN_CONFLICT : s->extra[v];
        s->num_conflicts++;
    }
    /* one of these voyages needs one more day */
    for (int v = 0; v < n; v++) {
        if (s->conflict[(size_t)hit * (size_t)n + (size_t)v] == NOT_IN_CONFLICT || s->extra[v] >= s->most[v])
            continue;
        s->extra[v]++;
        if (!seen_before(s, n) && vessels_fit(s, p, p->type[v], nullptr) && plan_search(s, p, berth, position))
            return true;
        s->extra[v]--;
        if (s->undecided)
            return false;
    }
    return false;
}

ScheduleStatus schedule_plan(Scheduler *s, const ScheduleProblem *p, int *vessel, int *berth, int *position)
{
    const int n = p->n, P = p->periods, D = p->day_length;
    s->stats.plan_calls++;
    if (n == 0)
        return SCHEDULE_FEASIBLE;
    RESERVE(s->base, s->cap_base, n);
    RESERVE(s->slack, s->cap_slack, n);
    RESERVE(s->reach, s->cap_reach, n);
    RESERVE(s->wait, s->cap_wait, n);
    RESERVE(s->trial, s->cap_trial, n);
    RESERVE(s->extra, s->cap_extra, n);
    RESERVE(s->most, s->cap_most, n);
    RESERVE(s->inside, s->cap_inside, n);
    RESERVE(s->type_start, s->cap_type_start, n);
    RESERVE(s->type_length, s->cap_type_length, n);
    RESERVE(s->type_index, s->cap_type_index, n);
    RESERVE(s->type_vessel, s->cap_type_vessel, n);
    for (int v = 0; v < n; v++) {
        /* days the vessel is busy without waiting, and the wait that fits in them */
        s->base[v] = p->time[v] / D + 1;
        if (s->base[v] > P)
            return SCHEDULE_INFEASIBLE;
        s->slack[v] = s->base[v] * D - p->time[v];
        s->extra[v] = 0;
    }
    for (int v = 0; v < n; v++) {
        /* the longest wait possible: a full carry plus the positions - 1 longest other loadings of the day */
        int64_t reach = p->max_carry;
        int last = -1; /* previous pick, in the order (longest first, lowest index first) */
        for (int k = 1; k < p->positions; k++) {
            int best = -1;
            for (int u = 0; u < n; u++) {
                if (u == v || p->day[u] != p->day[v])
                    continue;
                if (last >= 0 && (p->load[u] > p->load[last] || (p->load[u] == p->load[last] && u <= last)))
                    continue; /* picked already */
                if (best < 0 || p->load[u] > p->load[best])
                    best = u;
            }
            if (best < 0)
                break;
            reach += p->load[best];
            last = best;
        }
        s->reach[v] = reach > INT32_MAX / 2 ? INT32_MAX / 2 : (int)reach;
        /* extra days beyond which the limit can no longer bind, or that the week does not have */
        int most = 0;
        while ((int64_t)s->slack[v] + (int64_t)most * D < s->reach[v] && s->base[v] + most < P)
            most++;
        s->most[v] = (uint8_t)imin(most, NOT_IN_CONFLICT - 1);
    }
    for (int k = 0; k < p->num_types; k++)
        if (!vessels_fit(s, p, k, nullptr))
            return SCHEDULE_INFEASIBLE;
    s->num_seen = 0;
    s->num_seen_buckets = 0;
    s->num_conflicts = 0;
    s->nodes = 0;
    s->undecided = false;
    if (!plan_search(s, p, berth, position)) {
        s->stats.plan_limits += s->undecided;
        return s->undecided ? SCHEDULE_UNDECIDED : SCHEDULE_INFEASIBLE;
    }
    for (int v = 0; v < n; v++)
        if (s->extra[v] > 0) {
            s->stats.plan_waits++;
            break;
        }
    if (vessel)
        for (int k = 0; k < p->num_types; k++)
            if (!vessels_fit(s, p, k, vessel))
                die("schedule_plan: internal error (the vessels no longer fit)");
    return SCHEDULE_FEASIBLE;
}

/* ------------------------------------------------------------------------- */
/*  Direct check of an assignment                                            */
/* ------------------------------------------------------------------------- */

bool schedule_check(const ScheduleProblem *p, const int *vessel, const int *berth, const int *position)
{
    const int n = p->n, P = p->periods, B = p->berths, K = p->positions, D = p->day_length;
    bool ok = true;
    int64_t *load = XCALLOC_ARRAY(int64_t, P * B * K), *carry = XCALLOC_ARRAY(int64_t, P * B);
    bool *used = XCALLOC_ARRAY(bool, P * B * K);
    for (int v = 0; v < n && ok; v++) {
        ok = berth[v] >= 0 && berth[v] < B && position[v] >= 0 && position[v] < K && vessel[v] >= 0 &&
             vessel[v] < p->fleet[p->type[v]];
        if (!ok)
            break;
        const int slot = (p->day[v] * B + berth[v]) * K + position[v];
        ok = !used[slot];
        used[slot] = true;
        load[slot] = p->load[v];
    }
    /* carried waiting time: two rounds from idle berths reach the cyclic values, a third must agree */
    for (int b = 0; b < B && ok; b++) {
        int64_t w = 0;
        for (int round = 0; round < 3 && ok; round++)
            for (int d = 0; d < P && ok; d++) {
                if (round == 2)
                    ok = carry[d * B + b] == w;
                carry[d * B + b] = w;
                ok = ok && w <= p->max_carry;
                for (int k = 0; k < K; k++)
                    w += load[(d * B + b) * K + k];
                w = w > D ? w - D : 0;
            }
    }
    for (int v = 0; v < n && ok; v++) {
        int gap = P; /* days until the next departure of the vessel */
        for (int u = 0; u < n; u++)
            if (u != v && p->type[u] == p->type[v] && vessel[u] == vessel[v]) {
                const int g = (p->day[u] - p->day[v] + P) % P;
                gap = imin(gap, g);
            }
        int64_t wait = carry[p->day[v] * B + berth[v]];
        for (int k = 0; k < position[v]; k++)
            wait += load[(p->day[v] * B + berth[v]) * K + k];
        /* the vessel is busy on time / D + 1 days, and back in time after waiting */
        ok = gap >= p->time[v] / D + 1 && wait + p->time[v] <= (int64_t)gap * D;
    }
    free(load);
    free(carry);
    free(used);
    return ok;
}
