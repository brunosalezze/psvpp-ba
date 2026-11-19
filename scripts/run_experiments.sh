#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
# Runs the multi-start ALNS several times on a set of instances, as in the
# paper's computational experiments (Section 5.3: 15 runs per instance with the
# parameters of Table 3), and summarises best, average, deviation and CPU time.
#
# usage: scripts/run_experiments.sh [-r RUNS] [-j THREADS] [-o DIR] [-a "ALNS OPTIONS"] [INSTANCE...]
#
#   -r RUNS     runs per instance (default 15)
#   -j THREADS  runs executed in parallel (default: number of cores)
#   -o DIR      output directory (default results/<date-time>)
#   -a OPTIONS  extra options for `psvpp-ba alns`, e.g. -a "--iterations 50000 --uniform"
#   INSTANCE    instance names or files (default: the 12 non-clustered instances of Table 4)
#
# Output: DIR/<instance>.csv (one line per run) and DIR/summary.txt.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${PSVPP_BIN:-$ROOT/build/psvpp-ba}
RUNS=15
THREADS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
OUT=""
EXTRA=""

while [ $# -gt 0 ]; do
    case "$1" in
        -r) RUNS=$2; shift 2 ;;
        -j) THREADS=$2; shift 2 ;;
        -o) OUT=$2; shift 2 ;;
        -a) EXTRA=$2; shift 2 ;;
        -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
        *) break ;;
    esac
done

if [ $# -eq 0 ]; then
    set -- p10-1-1 p15-1-1 p21-1-1 p30-2-1 p41-2-1 p41-3-1 p50-3-1 p60-3-1 p69-4-1 p79-4-1 p79-5-1 p79-6-1
fi
[ -x "$BIN" ] || { echo "error: $BIN not found; run 'make' first" >&2; exit 1; }
OUT=${OUT:-$ROOT/results/$(date +%Y%m%d-%H%M%S)}
mkdir -p "$OUT"

printf '%-12s %5s %12s %12s %8s %9s\n' instance runs best average dev% cpu_s | tee "$OUT/summary.txt"
for inst in "$@"; do
    file=$inst
    [ -f "$file" ] || file=$ROOT/instances/$inst
    name=$(basename "$file")
    # shellcheck disable=SC2086 # EXTRA holds several options on purpose
    "$BIN" alns --csv -r "$RUNS" -j "$THREADS" $EXTRA "$file" > "$OUT/$name.csv"
    awk -F, -v name="$name" 'NR > 1 {
            n++; s += $3; q += $3 * $3; t += $6; if (n == 1 || $3 < b) b = $3; if ($5 != 1) inf++
        } END {
            a = s / n; sd = n > 1 ? sqrt((q - n * a * a) / (n - 1)) : 0
            printf "%-12s %5d %12.2f %12.2f %8.2f %9.1f%s\n", name, n, b, a, 100 * sd / a, t / n, inf ? "  (" inf " infeasible)" : ""
        }' "$OUT/$name.csv" | tee -a "$OUT/summary.txt"
done
echo "results written to $OUT"
