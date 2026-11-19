// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file solution.c
 * @brief Construction, copy and bookkeeping of heuristic solutions (see solution.h).
 */
#include "heuristic/solution.h"

#include "common/util.h"

#include <stdlib.h>
#include <string.h>

/** Violations up to this amount (seconds or load units) are rounding noise. */
enum { FEASIBILITY_TOLERANCE = 10 };

static size_t align8(size_t x) { return (x + 7u) & ~(size_t)7u; }

SolutionLayout sol_layout(const Instance *in)
{
    const size_t P = (size_t)in->num_periods, N = (size_t)in->num_installations;
    const size_t V = (size_t)in->num_vessels, B = (size_t)in->num_berths;
    SolutionLayout L = {0};
    size_t off = 0;
#define PLACE(field, bytes) (L.field = off, off = align8(off + (bytes)))
    PLACE(voyage, sizeof(Voyage) * V * P);
    PLACE(vessel, sizeof(Vessel) * V);
    PLACE(day, sizeof(BerthDay) * B * P);
    PLACE(berth, sizeof(BerthLog) * B);
    PLACE(node, sizeof(Node) * (size_t)in->num_nodes);
    PLACE(departure, sizeof(VoyageKey) * B * P * (size_t)in->slots_per_berth);
    PLACE(visit_node, sizeof(int) * P * N);
    PLACE(visit_count, sizeof(int8_t) * P * N);
    PLACE(pattern, sizeof(int8_t) * N);
    PLACE(active, sizeof(int8_t) * B * (size_t)in->num_types);
#undef PLACE
    L.size = off;
    return L;
}

/* Points the array fields of @p s into its memory block, following the layout @p L. */
static void bind_views(Solution *s, const SolutionLayout *L)
{
    unsigned char *b = s->block;
    s->voyage = (Voyage *)(b + L->voyage);
    s->vessel = (Vessel *)(b + L->vessel);
    s->day = (BerthDay *)(b + L->day);
    s->berth = (BerthLog *)(b + L->berth);
    s->node = (Node *)(b + L->node);
    s->departure = (VoyageKey *)(b + L->departure);
    s->visit_node = (int *)(b + L->visit_node);
    s->visit_count = (int8_t *)(b + L->visit_count);
    s->pattern = (int8_t *)(b + L->pattern);
    s->active = (int8_t *)(b + L->active);
}

Solution *sol_create(const Instance *inst)
{
    const SolutionLayout L = sol_layout(inst);
    Solution *s = XCALLOC_ARRAY(Solution, 1);
    s->inst = inst;
    s->block_size = L.size;
    s->block = xmalloc(L.size);
    bind_views(s, &L);
    sol_clear(s);
    return s;
}

void sol_free(Solution *s)
{
    if (s == nullptr)
        return;
    free(s->block);
    free(s);
}

void sol_copy(Solution *dst, const Solution *src)
{
    if (dst == src)
        return;
    dst->cost = src->cost;
    dst->infeasible = src->infeasible;
    dst->viol = src->viol;
    memcpy(dst->block, src->block, src->block_size);
}

/* Empty state of one vessel slot and its voyages. */
static void reset_vessel(Solution *s, int v)
{
    const Instance *in = s->inst;
    s->vessel[v] = (Vessel){.visits = 0, .first_period = -1, .excess_time = 0, .excess_load = 0, .cost = 0};
    for (int p = 0; p < in->num_periods; p++)
        *sol_voyage(s, v, p) = (Voyage){.size = 0, .queue_pos = -1, .next_period = -1, .prev_period = -1};
}

void sol_clear(Solution *s)
{
    const Instance *in = s->inst;
    memset(s->block, 0, s->block_size);
    s->cost = 0;
    s->infeasible = true;
    s->viol = (Violations){0};
    memset(s->pattern, -1, (size_t)in->num_installations);
    for (int h = in->num_visits; h < in->num_nodes; h++)
        s->node[h] = (Node){.next = h, .prev = h, .at = in->head_key[h - in->num_visits]};
    for (int v = 0; v < in->num_vessels; v++)
        reset_vessel(s, v);
}

void sol_update_cost(Solution *s, float w_load, float w_time, float w_berth)
{
    const Instance *in = s->inst;
    float total = 0;
    for (int b = 0; b < in->num_berths; b++) {
        const BerthLog *bl = &s->berth[b];
        total += (float)bl->excess_departures * w_berth + (float)bl->excess_time * w_time;
        for (int k = 0; k < in->num_types; k++)
            for (int l = 0; l < *sol_active(s, b, k); l++) {
                const Vessel *ve = sol_vessel(s, b, k, l);
                total += (float)ve->cost + (float)ve->excess_load * w_load + (float)ve->excess_time * w_time;
            }
    }
    s->cost = total;
}

double sol_objective(const Solution *s)
{
    const Instance *in = s->inst;
    double total = 0;
    for (int v = 0; v < in->num_vessels; v++) {
        if (s->vessel[v].visits == 0)
            continue;
        const VesselType *t = &in->type[in->head_key[v].type]; /* key of the day-0 head of slot v */
        total += t->fixed_cost;
        for (int p = 0; p < in->num_periods; p++)
            total += (double)t->var_cost * s->voyage[v * in->num_periods + p].dist;
    }
    return total / COST_SCALE;
}

void sol_check_feasibility(Solution *s)
{
    const Instance *in = s->inst;
    Violations n = {0};
    for (int b = 0; b < in->num_berths; b++) {
        const BerthLog *bl = &s->berth[b];
        if (bl->excess_departures > 0)
            n.berth++;
        if (bl->excess_time > FEASIBILITY_TOLERANCE)
            n.time++;
        for (int k = 0; k < in->num_types; k++)
            for (int l = 0; l < *sol_active(s, b, k); l++) {
                const Vessel *ve = sol_vessel(s, b, k, l);
                if (ve->excess_load > FEASIBILITY_TOLERANCE)
                    n.load++;
                if (ve->excess_time > FEASIBILITY_TOLERANCE)
                    n.time++;
            }
    }
    s->viol = n;
    s->infeasible = n.load > 0 || n.time > 0 || n.berth > 0;
}

/* Moves the last active slot of (berth, type) into slot @p slot, which must be
 * inactive, relinking its voyages to the head nodes of the new slot. */
static void move_last_vessel_into(Solution *s, int berth, int type, int slot)
{
    const Instance *in = s->inst;
    int8_t *active = sol_active(s, berth, type);
    const int last = --*active;
    if (slot == last)
        return;
    const int vd = inst_vessel(in, berth, type, slot), vs = inst_vessel(in, berth, type, last);
    s->vessel[vd] = s->vessel[vs];
    for (int p = 0; p < in->num_periods; p++) {
        Voyage *dst = sol_voyage(s, vd, p);
        const Voyage *src = sol_voyage(s, vs, p);
        *dst = *src;
        if (src->size > 0) {
            Node *nd = s->node;
            const int new_head = inst_head_node(in, p, vd), old_head = inst_head_node(in, p, vs);
            const int first = nd[old_head].next, lastn = nd[old_head].prev;
            nd[old_head].next = nd[old_head].prev = old_head;
            nd[new_head].next = first;
            nd[new_head].prev = lastn;
            nd[first].prev = new_head;
            nd[lastn].next = new_head;
            for (int u = first; u != new_head; u = nd[u].next)
                nd[u].at.slot = (int8_t)slot;
            sol_queue(s, berth, p)[dst->queue_pos].slot = (int8_t)slot;
        }
    }
    reset_vessel(s, vs);
}

void sol_compact(Solution *s)
{
    const Instance *in = s->inst;
    for (int b = 0; b < in->num_berths; b++)
        for (int k = 0; k < in->num_types; k++)
            for (int l = *sol_active(s, b, k) - 1; l >= 0; --l) {
                const int v = inst_vessel(in, b, k, l);
                if (s->vessel[v].visits == 0) {
                    move_last_vessel_into(s, b, k, l);
                    continue;
                }
                for (int p = in->num_periods - 1; p >= 0; --p) {
                    Voyage *r = sol_voyage(s, v, p);
                    if (r->size == 0) {
                        r->travel_time = 0;
                        r->time = 0;
                    }
                }
            }
}

void sol_fleet(const Solution *s, int *count)
{
    const Instance *in = s->inst;
    for (int k = 0; k < in->num_types; k++) {
        count[k] = 0;
        for (int b = 0; b < in->num_berths; b++)
            for (int l = 0; l < in->type[k].max_slots; l++)
                count[k] += s->vessel[inst_vessel(in, b, k, l)].visits > 0;
    }
}

/* ------------------------------------------------------------------------- */
/*  Printing                                                                 */
/* ------------------------------------------------------------------------- */

static double days(int seconds) { return (double)seconds / SECONDS_PER_DAY; }

void sol_print(const Solution *s, FILE *out)
{
    const Instance *in = s->inst;
    const int P = in->num_periods;
    int *fleet = XMALLOC_ARRAY(int, in->num_types);
    sol_fleet(s, fleet);

    double fixed = 0;
    for (int k = 0; k < in->num_types; k++)
        fixed += (double)fleet[k] * in->type[k].fixed_cost / COST_SCALE;
    const double total = sol_objective(s);
    fprintf(out, "objective %.6f (fleet %.6f + sailing %.6f)%s\n", total, fixed, total - fixed,
            s->infeasible ? "  [INFEASIBLE]" : "");
    fprintf(out, "fleet:");
    for (int k = 0; k < in->num_types; k++)
        fprintf(out, " type %d x %d%s", in->type[k].id, fleet[k], k + 1 < in->num_types ? "," : "\n");
    if (s->infeasible)
        fprintf(out, "violations: %d load, %d time, %d berth\n", s->viol.load, s->viol.time, s->viol.berth);

    fprintf(out, "voyages (day: berth vessel type#slot | load m2 | wait+load+sail = total days | route):\n");
    int count = 0;
    for (int p = 0; p < P; p++)
        for (int b = 0; b < in->num_berths; b++) {
            const int nq = sol_day(s, b, p)->num_departures;
            const VoyageKey *q = sol_queue(s, b, p);
            for (int j = 0; j < nq; j++) {
                const Voyage *r = sol_voyage_of(s, q[j]);
                fprintf(out, "  day %d: berth %d vessel %d#%d | %7.2f | %.2f+%.2f+%.2f = %.2f | berth", p + 1,
                        in->berth[b].id, in->type[q[j].type].id, q[j].slot + 1, r->load / 1e3, days(r->wait_time),
                        days(r->load_time), days(r->travel_time), days(r->time));
                const int head = inst_head_node(in, p, sol_vessel_of(s, q[j]));
                for (int u = s->node[head].next; u != head; u = s->node[u].next)
                    fprintf(out, " -> %d", in->inst[in->visit_inst[u]].id);
                fprintf(out, " -> berth\n");
                count++;
            }
        }
    fprintf(out, "%d voyages\n", count);

    fprintf(out, "berth schedules (day: carried wait + loading = occupation, days):\n");
    for (int b = 0; b < in->num_berths; b++) {
        fprintf(out, "  berth %d:", in->berth[b].id);
        for (int p = 0; p < P; p++) {
            const BerthDay *d = sol_day(s, b, p);
            fprintf(out, " [%d] %.2f+%.2f", p + 1, days(d->carry_time - d->sink), days(d->load_time));
        }
        fputc('\n', out);
    }
    fprintf(out, "delivery patterns (installation: days):\n ");
    for (int i = 0; i < in->num_installations; i++) {
        fprintf(out, " %d:", in->inst[i].id);
        for (int p = 0; p < P; p++)
            fputc(s->visit_count[p * in->num_installations + i] ? (P <= 9 ? '1' + p : 'x') : '.', out);
        if ((i + 1) % 8 == 0 && i + 1 < in->num_installations)
            fprintf(out, "\n ");
    }
    fputc('\n', out);
    free(fleet);
}
