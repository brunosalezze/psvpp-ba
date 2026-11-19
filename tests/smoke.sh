#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
# End-to-end smoke tests of a built psvpp-ba binary (no Gurobi needed).
#
#   ./tests/smoke.sh build/psvpp-ba
#
# Runs the Gurobi-free commands on the small instances, checks what they print,
# has scripts/check_solution.py verify the solutions the heuristic saves, and
# checks that the exact commands of a build without Gurobi report that instead
# of crashing. Prints a line per step and "smoke tests OK" at the end; exits
# non-zero on the first unexpected result.
set -euo pipefail

usage() {
    printf 'usage: %s <path-to-psvpp-ba>\n' "$0" >&2
    printf 'example: make test    # or ./tests/smoke.sh build/psvpp-ba\n' >&2
}

if [ "$#" -ne 1 ]; then
    usage
    exit 2
fi

# the paths below are relative to the repository, wherever the script is called from
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

# accept the binary as a path of the caller or of the repository (make test runs from the root)
bin=$1
if [ ! -f "$bin" ] || [ ! -x "$bin" ]; then
    if [ -f "$root_dir/$bin" ] && [ -x "$root_dir/$bin" ]; then
        bin=$root_dir/$bin
    else
        printf 'error: %s is not an executable file\n' "$1" >&2
        usage
        exit 2
    fi
elif [ "${bin:0:1}" != "/" ]; then
    bin=$PWD/$bin
fi
CDPATH= cd -- "$root_dir"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/psvpp-smoke.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

instances=(p10-1-1 p15-1-1 p30-2-1 p30-2-4)
# a quick but non-trivial ALNS run (start-iterations defaults to 100)
alns_opts=(--starts 2 --iterations 2000 --seed 1)

step() { printf '== %s\n' "$*"; }

# fail <message>: report what went wrong and stop; the first difference is the one to look at
fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

dump() { sed 's/^/    | /' "$1" >&2 || true; }

# capture <file in $tmp> <command...>: run it, keep stdout and stderr, ignore the status
capture() {
    local cap=$tmp/$1
    shift
    "$@" >"$cap" 2>&1 || true
    return 0
}

# expect_ok <what> <file in $tmp> <command...>: the command must succeed
expect_ok() {
    local what=$1 cap=$tmp/$2 status=0
    shift 2
    "$@" >"$cap" 2>&1 || status=$?
    if [ "$status" -ne 0 ]; then
        dump "$cap"
        fail "$what: '$*' exited with status $status"
    fi
    return 0
}

# expect_fail <what> <file in $tmp> <pattern> <command...>: the command must fail mentioning pattern
expect_fail() {
    local what=$1 cap=$tmp/$2 pattern=$3
    shift 3
    if "$@" >"$cap" 2>&1; then
        fail "$what: '$*' was expected to fail but succeeded"
    fi
    grep -qE -- "$pattern" "$cap" || {
        dump "$cap"
        fail "$what: expected the output to match '$pattern'"
    }
    return 0
}

# grep_for <what> <file in $tmp> <pattern>: the captured output must match pattern
grep_for() {
    local what=$1 cap=$tmp/$2 pattern=$3
    grep -qE -- "$pattern" "$cap" || {
        dump "$cap"
        fail "$what: expected the output to match '$pattern'"
    }
    return 0
}

# grep_not <what> <file in $tmp> <pattern>: the captured output must not match pattern
grep_not() {
    local what=$1 cap=$tmp/$2 pattern=$3
    if grep -qE -- "$pattern" "$cap"; then
        dump "$cap"
        fail "$what: expected no output matching '$pattern'"
    fi
    return 0
}

# strip_cpu <file in $tmp> <file in $tmp>: copy removing the CPU time of the result lines
strip_cpu() { sed 's/, cpu .*/, cpu/' "$tmp/$1" >"$tmp/$2"; }

# ---------------------------------------------------------------------------
step "version and usage"
expect_ok "--version" version "$bin" --version
grep_for "--version" version '^psvpp-ba [0-9]'
expect_ok "--help" help "$bin" --help
for command in alns routes bc complete info; do
    grep_for "--help lists $command" help "(^| )$command( |\$)"
done
expect_fail "no arguments" no-args 'usage: psvpp-ba' "$bin"
expect_fail "unknown command" bad-cmd 'unknown command' "$bin" nosuchcommand instances/p10-1-1
expect_fail "missing instance" no-instance 'exactly one instance file' "$bin" alns
expect_fail "unreadable instance" bad-instance 'cannot open instance file' "$bin" alns "$tmp/no-such-instance"
expect_fail "impossible iterations" bad-iters 'smaller than starts' \
    "$bin" alns --starts 2 --iterations 100 instances/p10-1-1

step "info"
for inst in "${instances[@]}"; do
    expect_ok "info $inst" "info-$inst" "$bin" info "instances/$inst"
    grep_for "info $inst" "info-$inst" "^instance instances/$inst: [1-9][0-9]* installations"
done

# ---------------------------------------------------------------------------
step "alns on ${instances[*]}"
for inst in "${instances[@]}"; do
    sol=$tmp/solution-$inst
    expect_ok "alns $inst" "alns-$inst" "$bin" alns "${alns_opts[@]}" --output "$sol" "instances/$inst"
    # one result line per run, feasible, with an objective and a fleet
    grep_for "alns $inst (result line)" "alns-$inst" 'run 1 \(seed 1\): objective [0-9]+\.[0-9]+, fleet \[[0-9 ]+\]'
    grep_not "alns $inst (feasibility)" "alns-$inst" 'INFEASIBLE'
    # the best solution is printed and looks like a plan
    grep_for "alns $inst (best solution)" "alns-$inst" '^best solution:$'
    grep_for "alns $inst (objective)" "alns-$inst" \
        '^objective [0-9]+\.[0-9]+ \(fleet [0-9]+\.[0-9]+ \+ sailing [0-9]+\.[0-9]+\)$'
    grep_for "alns $inst (voyages)" "alns-$inst" '^  day [0-9]+: berth [0-9]+ vessel [0-9]+#[0-9]+ \|'
    grep_for "alns $inst (voyage count)" "alns-$inst" '^[0-9]+ voyages$'
    grep_for "alns $inst (berth schedule)" "alns-$inst" '^berth schedules'
    # --output keeps the same solution as stdout
    [ -s "$sol" ] || fail "alns $inst: --output wrote nothing to $sol"
    printed=$(sed -n 's/^objective \([0-9.]*\) (fleet .*/\1/p' "$tmp/alns-$inst" | tail -n 1)
    saved=$(sed -n 's/^objective \([0-9.]*\) (fleet .*/\1/p' "$sol" | head -n 1)
    [ -n "$printed" ] && [ "$printed" = "$saved" ] ||
        fail "alns $inst: the objective in $sol ($saved) is not the printed one ($printed)"
done

step "alns --quiet prints only the result lines"
capture quiet "$bin" alns "${alns_opts[@]}" --quiet instances/p10-1-1
grep_for "alns --quiet" quiet '^run 1 \(seed 1\): objective '
grep_not "alns --quiet" quiet 'best solution'
grep_not "alns --quiet" quiet '^instance instances/'

step "alns --csv with several runs (in parallel when OpenMP is available)"
expect_ok "alns --csv" csv "$bin" alns --runs 3 --threads 3 "${alns_opts[@]}" --csv instances/p30-2-1
grep_for "alns --csv (header)" csv '^instance,seed,objective,penalised,feasible,cpu_s,fleet_type'
rows=$(grep -c '^instances/p30-2-1,' "$tmp/csv" || true)
[ "$rows" -eq 3 ] || fail "alns --csv: expected 3 result rows, found $rows"
# the feasibility column (5th) must be 1 for every run; a build without OpenMP also
# prints a warning on stderr, which the capture mixes with the rows, so look at the rows only
infeasible=$(awk -F, '$1 ~ /^instances\// && $5 != 1 { n++ } END { print n + 0 }' "$tmp/csv")
[ "$infeasible" -eq 0 ] || fail "alns --csv: $infeasible of 3 runs were infeasible"

step "alns --verbose --uniform (uniform operator choice, operator statistics)"
expect_ok "alns --verbose" verbose "$bin" alns --starts 1 --iterations 500 --ls-interval 25 --verbose --uniform \
    instances/p10-1-1
grep_for "alns --verbose" verbose 'uniform operator selection'
grep_for "alns --verbose (operators)" verbose 'operator +uses +final weight'

step "the same seed gives the same result (strict floating point)"
for twice in 1 2; do
    capture "run-$twice" "$bin" alns --starts 2 --iterations 2000 --seed 7 --quiet instances/p15-1-1
    strip_cpu "run-$twice" "run-$twice.norm"
done
cmp "$tmp/run-1.norm" "$tmp/run-2.norm" || fail "alns: two runs with seed 7 printed different results"

step "a different seed explores a different solution"
capture seed-8 "$bin" alns --starts 2 --iterations 2000 --seed 8 --quiet instances/p15-1-1
strip_cpu seed-8 seed-8.norm
if cmp -s "$tmp/run-1.norm" "$tmp/seed-8.norm"; then
    fail "alns: seeds 7 and 8 gave exactly the same run"
fi

# ---------------------------------------------------------------------------
step "scripts/check_solution.py verifies the saved solutions"
if ! command -v python3 >/dev/null 2>&1; then
    printf '   (python3 not found: skipping the validator)\n'
else
    for inst in "${instances[@]}"; do
        # the heuristic does not restrict a voyage to one cluster (see the note in
        # src/common/instance.h: only the voyage enumeration of the exact models does),
        # so a plan of `alns` on a clustered instance is not a plan of the model.
        clusters=$(sed -n 's/.* \([0-9][0-9]*\) cluster.*/\1/p' "$tmp/info-$inst")
        if [ "${clusters:-1}" -gt 1 ]; then
            printf '   (%s: %s clusters, the heuristic ignores them: not validated)\n' "$inst" "$clusters"
            continue
        fi
        expect_ok "check_solution.py $inst" "check-$inst" python3 scripts/check_solution.py \
            "instances/$inst" "$tmp/solution-$inst"
        grep_for "check_solution.py $inst" "check-$inst" 'solution is valid$'
        grep_not "check_solution.py $inst" "check-$inst" 'ERROR:'
    done
fi

# ---------------------------------------------------------------------------
step "routes (voyage enumeration, needs no Gurobi)"
for inst in p10-1-1 p15-1-1 p30-2-1; do
    expect_ok "routes $inst" "routes-$inst" "$bin" routes "instances/$inst"
    grep_for "routes $inst" "routes-$inst" '^voyages: [1-9][0-9]* \('
    grep_for "routes $inst (histogram)" "routes-$inst" '^voyages per number of installations: 1:[0-9]+'
    grep_for "routes $inst (distances)" "routes-$inst" '^sum of voyage distances: [0-9]+$'
done
# --max-route-size really limits the voyages that are enumerated
expect_ok "routes --max-route-size" routes-limited "$bin" routes --max-route-size 4 instances/p15-1-1
grep_for "routes --max-route-size" routes-limited \
    '^voyages per number of installations: 1:[0-9]+ 2:[0-9]+ 3:[0-9]+ 4:[0-9]+$'

# ---------------------------------------------------------------------------
if LC_ALL=C grep -aqF 'no Gurobi support' "$bin"; then
    step "bc and complete say that this build has no Gurobi"
    for command in bc complete; do
        expect_fail "$command without Gurobi" "no-gurobi-$command" 'no Gurobi support' \
            "$bin" "$command" --time-limit 1 instances/p10-1-1
    done
else
    step "this build has Gurobi: the exact commands are not exercised here"
fi

# ---------------------------------------------------------------------------
printf 'smoke tests OK\n'
