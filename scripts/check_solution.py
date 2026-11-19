#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
"""Checks a solution of the exact methods against the instance file, independently of the C program.

usage: scripts/check_solution.py [--max-route-size N] INSTANCE SOLUTION

SOLUTION is the file written by `psvpp-ba bc -o FILE`, `psvpp-ba complete -o FILE` or
`psvpp-ba alns -o FILE` (the "best solution" printed by these commands can be saved to a file as
well). In a solution of the heuristic a vessel is based at a berth and named by its type and its
slot there. Everything is recomputed from the instance data (docs/instance-format.md) and the rules
of the complete model of the paper:

  * every installation is visited on the days of one of its delivery patterns;
  * a voyage stays in one cluster, visits at most 8 installations (--max-route-size, as the option
    of the program) and its load fits in its vessel type (it may fill it exactly);
  * a berth has at most two departures per day, loaded one after the other; the loading that does not
    fit in the day is carried to the next day as waiting time, at most half a day, the week being cyclic;
  * a vessel is back from a voyage (waiting for the berth included) before its next departure;
  * the fleet and the cost are the ones announced, and the fleet respects the limit of the method:
    vessels per type for the exact models, vessels per type and berth for the heuristic.

Times are computed in double precision from the input, whereas the program uses integer seconds with
the rounding of the original code; the smallest margins are printed so that a difference of a few
seconds can be told from a real violation. Exit status 0 if the solution passes every check.
"""
import math
import re
import sys

DAY = 86400.0
TOLERANCE = 1.5  # days a berth may be occupied per day
POSITIONS = 2    # departures per berth and day
SLACK = 60.0     # seconds: differences below this are rounding
LOAD_SLACK = 0.005  # m2: the program rounds the load of a visit to 0.001 m2


def read_instance(path):
    tok = open(path).read().split()
    pos = 0

    def take(n):
        nonlocal pos
        out = [float(x) for x in tok[pos:pos + n]]
        pos += n
        return out

    T, B, N, P, C = (int(x) for x in take(5))
    types = {}
    for _ in range(T):
        f = take(9)
        types[int(f[0])] = dict(slots=int(f[1]), max_fleet=int(f[2]), fixed=f[3], var=f[4], capacity=f[5], knots=f[8])
    berths = [dict(x=f[1], y=f[2], setup=f[3] * DAY, rate=f[4] * DAY) for f in (take(5) for _ in range(B))]
    inst = {}
    for _ in range(N):
        f = take(11)
        inst[int(f[0])] = dict(cluster=int(f[1]), x=f[2], y=f[3], demand=f[4], service=f[6] * DAY,
                               unload=f[7] * DAY / 1000.0, freq=int(f[8]))
    patterns = {}
    for _ in range(int(take(1)[0])):
        freq, count = (int(x) for x in take(2))
        for code in take(count):
            days = frozenset(d for d in range(P) if (int(code) >> (P - 1 - d)) & 1)
            patterns.setdefault(freq, set()).add(days)
    return dict(T=T, B=B, N=N, P=P, C=C, types=types, berths=berths, inst=inst, patterns=patterns)


def read_solution(path):
    text = open(path).read()
    m = re.search(r"objective ([\d.]+) \(fleet ([\d.]+) \+ sailing ([\d.]+)\)", text)
    if not m:
        sys.exit(f"{path}: no 'objective' line")
    objective = float(m.group(1))
    fleet = {int(a): int(b) for a, b in re.findall(r"type (\d+) x (\d+)", text)}
    voyages = []
    for m in re.finditer(r"day (\d+): berth (\d+)/(\d+) vessel (\d+)#(\d+) \|\s*([\d.]+) \| ([\d.]+) \+ ([\d.]+) \+ ([\d.]+) \| berth((?: -> \d+)+) -> berth", text):
        stops = [int(x) for x in re.findall(r"-> (\d+)", m.group(10))]
        voyages.append(dict(day=int(m.group(1)) - 1, berth=int(m.group(2)), position=int(m.group(3)),
                            type=int(m.group(4)), vessel=int(m.group(5)), load=float(m.group(6)),
                            wait=float(m.group(7)) * DAY, stops=stops))
    if voyages:
        return objective, fleet, voyages, False
    # solution of the heuristic: "day D: berth B vessel TYPE#SLOT | load | wait+load+sail = total | route";
    # the vessel is (berth, type, slot) and the voyages of a berth and day load in the order of their waiting times
    for m in re.finditer(r"day (\d+): berth (\d+) vessel (\d+)#(\d+) \|\s*([\d.]+) \| ([\d.]+)\+([\d.]+)\+([\d.]+) = [\d.]+ \| berth((?: -> \d+)+) -> berth", text):
        stops = [int(x) for x in re.findall(r"-> (\d+)", m.group(9))]
        voyages.append(dict(day=int(m.group(1)) - 1, berth=int(m.group(2)), type=int(m.group(3)),
                            vessel=1000 * int(m.group(2)) + int(m.group(4)), slot=int(m.group(4)), load=float(m.group(5)),
                            wait=float(m.group(6)) * DAY, stops=stops))
    order = {}
    for v in sorted(voyages, key=lambda v: v["wait"]):
        v["position"] = order[(v["day"], v["berth"])] = order.get((v["day"], v["berth"]), 0) + 1
    return objective, fleet, voyages, True


def main():
    max_stops = 8
    if len(sys.argv) > 2 and sys.argv[1] == "--max-route-size":
        max_stops = int(sys.argv[2])
        del sys.argv[1:3]
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    I = read_instance(sys.argv[1])
    objective, fleet, voyages, heuristic = read_solution(sys.argv[2])
    notes = []
    P, B = I["P"], I["B"]
    errors = []
    margins = {}

    def margin(name, value):
        margins[name] = min(margins.get(name, math.inf), value)

    base = I["berths"][0]
    knots = next(iter(I["types"].values()))["knots"]
    metre_seconds = 3.6 / (1.852 * knots)

    def dist(a, b):  # 1000 x Euclidean distance, rounded as in the program
        return math.floor(1000.0 * math.hypot(a["x"] - b["x"], a["y"] - b["y"]) + 0.5)

    # ---- voyages: cluster, load, times, cost
    cost = sum(I["types"][k]["fixed"] * n for k, n in fleet.items())
    visits = {i: set() for i in I["inst"]}
    for v in voyages:
        stops = [I["inst"][i] for i in v["stops"]]
        if I["C"] > 1 and len({s["cluster"] for s in stops}) > 1:
            errors.append(f"day {v['day'] + 1}: voyage {v['stops']} mixes clusters")
        if len(stops) > max_stops:
            errors.append(f"day {v['day'] + 1}: voyage {v['stops']} visits more than {max_stops} installations")
        load = sum(s["demand"] / s["freq"] for s in stops)
        capacity = I["types"][v["type"]]["capacity"]
        if load > capacity + LOAD_SLACK:
            errors.append(f"day {v['day'] + 1}: load {load:.2f} of {v['stops']} exceeds the capacity {capacity:g}")
        if abs(load - v["load"]) > 0.01:
            errors.append(f"day {v['day'] + 1}: load {load:.2f} differs from the announced {v['load']:.2f}")
        points = [base] + stops + [base]
        length = sum(dist(points[k], points[k + 1]) for k in range(len(points) - 1))
        v["loading"] = base["setup"] + sum(s["demand"] / s["freq"] * base["rate"] for s in stops)
        v["time"] = v["loading"] + length * metre_seconds + sum(
            s["service"] + s["demand"] / s["freq"] * s["unload"] for s in stops)
        cost += I["types"][v["type"]]["var"] * length / 1000.0  # the variable cost is per unit of coordinates
        for i in v["stops"]:
            if v["day"] in visits[i]:
                errors.append(f"installation {i} is visited twice on day {v['day'] + 1}")
            visits[i].add(v["day"])

    # ---- delivery patterns
    for i, days in visits.items():
        freq = I["inst"][i]["freq"]
        if frozenset(days) not in I["patterns"].get(freq, set()):
            errors.append(f"installation {i} (frequency {freq}) is visited on days "
                          f"{sorted(d + 1 for d in days)}, which is not one of its patterns")

    # ---- berths: positions, carried waiting time (cyclic), waiting time of every voyage
    slots = {}
    for v in voyages:
        key = (v["day"], v["berth"], v["position"])
        if key in slots or not 1 <= v["position"] <= POSITIONS:
            errors.append(f"day {v['day'] + 1}: position {v['position']} of berth {v['berth']} is used twice or invalid")
        slots[key] = v
    if len({v["berth"] for v in voyages}) > B:
        errors.append(f"more than {B} berths are used")
    for b in {v["berth"] for v in voyages}:
        carry = 0.0
        for sweep in range(4):  # the carried time settles after two turns of the week
            for d in range(P):
                if sweep == 3:
                    margin("berth carry below half a day", (TOLERANCE - 1.0) * DAY - carry)
                    ahead = carry
                    for k in range(1, POSITIONS + 1):
                        if (d, b, k) in slots:
                            v = slots[(d, b, k)]
                            if abs(ahead - v["wait"]) > 0.006 * DAY:
                                errors.append(f"day {d + 1}, berth {b}: waiting time {ahead / DAY:.3f} d differs from "
                                              f"the announced {v['wait'] / DAY:.2f} d")
                            v["true_wait"] = ahead
                            ahead += v["loading"]
                occupied = carry + sum(slots[(d, b, k)]["loading"] for k in range(1, POSITIONS + 1) if (d, b, k) in slots)
                carry = max(0.0, occupied - DAY)
        if carry > (TOLERANCE - 1.0) * DAY + SLACK:
            errors.append(f"berth {b} carries {carry / DAY:.3f} d of waiting time into a day")

    # ---- vessels: back before the next departure
    by_vessel = {}
    for v in voyages:
        by_vessel.setdefault((v["type"], v["vessel"]), []).append(v)
    for (k, n), trips in by_vessel.items():
        trips.sort(key=lambda v: v["day"])
        for a, v in enumerate(trips):
            nxt = trips[(a + 1) % len(trips)]
            gap = (nxt["day"] - v["day"] - 1) % P + 1 if len(trips) > 1 else P
            if len(trips) > 1 and nxt["day"] == v["day"]:
                errors.append(f"vessel {k}#{n} departs twice on day {v['day'] + 1}")
                continue
            margin("vessel back before its next departure", gap * DAY - (v.get("true_wait", v["wait"]) + v["time"]))
    for k in I["types"]:
        used = len({n for (t, n) in by_vessel if t == k})
        allowed = I["types"][k]["max_fleet"]
        if used > fleet.get(k, 0) or (not heuristic and fleet.get(k, 0) > allowed):
            errors.append(f"type {k}: {used} vessels used, {fleet.get(k, 0)} announced, {allowed} allowed")
        if heuristic:
            slots = I["types"][k]["slots"]
            if any(v["slot"] > slots for v in voyages if v["type"] == k):
                errors.append(f"type {k}: a berth bases more than its {slots} vessels of this type")
            if fleet.get(k, 0) > allowed:
                notes.append(f"type {k}: {fleet.get(k, 0)} vessels, more than the {allowed} of the exact models (max_fleet)")

    for name, value in margins.items():
        if value < -SLACK:
            errors.append(f"{name}: violated by {-value:.0f} s")
    if abs(cost - objective) > 1e-3:
        errors.append(f"cost {cost:.6f} differs from the announced objective {objective:.6f}")

    print(f"{sys.argv[2]}: {len(voyages)} voyages, {sum(len(v['stops']) for v in voyages)} visits, "
          f"fleet {dict(sorted(fleet.items()))}, cost {cost:.6f}")
    for name, value in margins.items():
        print(f"  smallest margin, {name}: {value:.0f} s")
    for n in notes:
        print("  note:", n)
    for e in errors:
        print("  ERROR:", e)
    print("  solution is valid" if not errors else f"  {len(errors)} error(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
