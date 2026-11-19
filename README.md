# psvpp-ba

Exact and heuristic algorithms for the **periodic supply vessel planning problem
with berth allocation (PSVPP-BA)**: choose the fleet of offshore supply vessels
(how many of each type), their weekly voyages to the offshore oil and gas
installations, and the berth schedule at the supply base, at minimum fixed and
sailing cost. The planning week repeats indefinitely, installations must be
visited according to their delivery patterns, and each berth loads at most two
vessels per day within a 1.5-day tolerance.

This is the implementation accompanying the paper

> B. S. Vieira, G. M. Ribeiro, L. Bahiense, R. Cruz, A. B. Mendes, G. Laporte.
> *Exact and heuristic algorithms for the fleet composition and periodic routing
> problem of offshore supply vessels with berth allocation decisions.*
> European Journal of Operational Research 295(3):908–923, 2021.
> [doi:10.1016/j.ejor.2021.03.022](https://doi.org/10.1016/j.ejor.2021.03.022)

It provides

* `alns` — the multi-start **adaptive large neighbourhood search** with spaced
  local searches (paper, Section 4); no external dependency;
* `bc` — the **branch-and-cut** on the reduced model (Sections 2.2 and 3);
* `complete` — the **complete model** (Section 2.1);
* `routes` — the enumeration of the voyages used by the exact models.

The exact methods use [Gurobi](https://www.gurobi.com) (free academic licences
are available); everything else builds with a C compiler alone.

## Building

Requirements: a C compiler supporting C23 (tested with GCC 13.2, 14.2, 14.3
and 15.2, Clang 23.1.5 and Intel oneAPI `icx` 2025.0, which build the code
without warnings and give identical results for the same seed, on Intel and
AMD processors), GNU make, and optionally OpenMP (parallel runs) and Gurobi
(exact methods; tested with Gurobi 13.0).

```sh
make            # optimised build: build/psvpp-ba
make test       # unit tests and end-to-end smoke tests
```

The build reports what it found, e.g.
`built build/psvpp-ba (C standard: c23, Gurobi: gurobi130, OpenMP: yes)`.

| variable | effect |
|---|---|
| `GUROBI_HOME=/path/to/gurobi/linux64` | where Gurobi is installed (the Gurobi installer sets it); when it is unset or names a removed version (with a warning), the newest `/opt/gurobi*/linux64` is used |
| `GUROBI=0` | build without Gurobi: `alns` and `routes` only |
| `OPENMP=0` | build without OpenMP: runs are executed one after the other |
| `CC=clang` | choose the compiler |
| `LDFLAGS=...` | extra linker flags (see the Clang note below) |
| `CFLAGS="-O3 -march=native"` | change optimisation flags |

`make debug` builds `build-debug/psvpp-ba` with AddressSanitizer and
UndefinedBehaviorSanitizer. `make help` lists all targets.

**Linux.** Install `gcc` (or `clang`) and `make` from the distribution
(`sudo apt install build-essential` on Debian/Ubuntu; GCC 13+ is in Debian 13,
Ubuntu 24.04 and later). OpenMP comes with GCC.

**Clang.** `make CC=clang` works with a Clang installed by the system package
manager (on Debian/Ubuntu also install `libomp-dev` for OpenMP). With a
standalone LLVM release unpacked in some directory, tell the linker where its
OpenMP runtime is:
`make CC=$LLVM/bin/clang LDFLAGS=-Wl,-rpath,$LLVM/lib/x86_64-unknown-linux-gnu`.

**macOS.** Apple Clang lacks OpenMP: either `make OPENMP=0`, or install GCC with
Homebrew (`brew install gcc`) and build with `make CC=gcc-15` (use the version
installed). Gurobi is found in
`/Library/gurobi*/macos_universal2` or through `GUROBI_HOME`.

**Windows.** Use WSL (Windows Subsystem for Linux) and follow the Linux steps.

**Gurobi licence.** The program reads the licence like any Gurobi application
(`~/gurobi.lic` or `GRB_LICENSE_FILE`). Run `gurobi_cl --license` to check it.

## Running

```sh
build/psvpp-ba <command> [options] <instance-file>
build/psvpp-ba --help                       # all options with their defaults
```

### Heuristic

```sh
# one run with the tuned parameters of the paper (alpha=50, beta=100, gamma=300,000)
build/psvpp-ba alns instances/p10-1-1

# 15 runs with seeds 1..15, 8 at a time, results as CSV
build/psvpp-ba alns --runs 15 --threads 8 --csv instances/p41-3-1

# a quick run: fewer starts and iterations
build/psvpp-ba alns --starts 10 --iterations 20000 instances/p79-4-1
```

The output gives, per run, the objective (fixed + sailing cost), the fleet
(vessels of each type) and the CPU time; with several runs, a summary (best,
average, standard deviation in %, average CPU time); and the best solution:
every voyage with its day, berth, vessel, load, waiting/loading/sailing times
and route, the berth schedules and the delivery pattern of each installation.
`--output FILE` saves the best solution; `--verbose` also prints how often each
operator was used, with its final adaptive weight, and, when runs do not
execute in parallel, the result of every start; `--quiet` prints only the
result lines. `--seed N` sets the seed of the first run (run *r* uses
`N + r - 1`); a seed fully determines a run, also when runs are executed in
parallel.

### Exact methods

```sh
build/psvpp-ba routes instances/p15-1-1                     # voyage set only
build/psvpp-ba bc instances/p10-1-1                         # branch-and-cut
build/psvpp-ba bc --time-limit 3600 --solver-log instances/p15-1-2
build/psvpp-ba complete --time-limit 3600 instances/p10-1-1 # complete model
```

The result line reports the best solution (UB), the bound (LB), the gap and the
solver, wall-clock and CPU times, followed by the solution. A solution saved
with `--output FILE` can be verified from the instance data alone with
`scripts/check_solution.py INSTANCE FILE`, a small Python program that shares
no code with the solver; it reads the solutions of `alns --output` as well.
`bc` rewrites that file whenever it accepts a better plan, and prints the bound
of the root LP relaxation at once and a progress line every minute, so a run
killed by a memory or time limit still leaves its best plan and bound.
`--solver-threads`
sets the threads of the main Gurobi model (default: all cores),
`--subproblem-threads` those of the small models that the branch-and-cut
falls back on when its scheduling algorithms cannot settle an incumbent
(default 1; see [docs/algorithms.md](docs/algorithms.md)), and
`--max-route-size` the largest number of installations per voyage (default 8,
as in the paper; the heuristic follows the same option). `--gap` sets the gap at
which the exact methods stop (default 0.01 %, Gurobi's; the original branch-and-cut
stopped at 0.1 %). `--min-usage 60` reproduces a filter of
the original implementation, which dropped the voyages loading less than 60 % of
the vessel on non-clustered instances: much smaller models, but no longer an
exact method (the default, 0, keeps every voyage). The number of voyages grows very
quickly with the number of installations of a cluster (about 0.8 million for
41 non-clustered installations and 29 million for 79), so the exact methods
are practical on the clustered instances and the smaller non-clustered ones.

### Reproducing the experiments

```sh
scripts/run_experiments.sh                     # 15 runs on the 12 instances of Table 4
scripts/run_experiments.sh -r 15 p10-1-1 p21-1-1
scripts/run_experiments.sh -a "--iterations 50000" p41-3-1
```

The script writes one CSV per instance and a `summary.txt` under `results/`.

## Instances

The 30 instances of the paper are in [`instances/`](instances) (licensed
under CC BY 4.0; [instances/README.md](instances/README.md) lists their origin);
the file format is described in [docs/instance-format.md](docs/instance-format.md).
`p<N>-<B>-<C>` has `N` installations, `B` berths and `C` clusters
(`C = 1`: non-clustered).

## Parameters

Defaults are the tuned values of the paper (Table 3).

| option | symbol | default | meaning |
|---|---|---|---|
| `--starts` | $\alpha$ | 50 | number of initial solutions |
| `--start-iterations` | $\beta$ | 100 | ALNS iterations per start |
| `--iterations` | $\gamma$ | 300,000 | total ALNS iterations |
| `--ls-interval` | $\delta$ | 50 | iterations between local searches |
| `--ls-max-rounds` | | 80 | improvement steps per local search, at most |
| `--segment` | $\rho$ | 50 | iterations between weight updates |
| `--t0`, `--tf` | $\theta'_0, \theta'_F$ | 1/50, 1/500 | initial/final temperature relative to the cost |
| `--removal-min`, `--removal-max` | $n^-_1, n^-_2$ | 0.10, 0.25 | fraction of installations removed |
| `--determinism` | $p$ | 4 | determinism of the guided removals |
| `--reaction` | $\tau_A$ | 0.05 | reaction factor of the adaptive layer |
| `--scores` | $\sigma_1,\sigma_2,\sigma_3$ | 5,0,10 | scores for new best / improving / accepted |
| `--penalty-increase`, `--penalty-decrease` | | 1.15, 0.90 | penalty weight update factors |
| `--uniform` | | off | choose operators uniformly instead of adaptively |

## Citation

If you use this code, please cite the paper
([CITATION.cff](CITATION.cff) holds the same reference in machine-readable form):

```bibtex
@article{Vieira2021,
  author  = {Vieira, Bruno S. and Ribeiro, Glaydston M. and Bahiense, Laura and
             Cruz, Roberto and Mendes, Andr{\'e} B. and Laporte, Gilbert},
  title   = {Exact and heuristic algorithms for the fleet composition and periodic
             routing problem of offshore supply vessels with berth allocation decisions},
  journal = {European Journal of Operational Research},
  volume  = {295},
  number  = {3},
  pages   = {908--923},
  year    = {2021},
  doi     = {10.1016/j.ejor.2021.03.022}
}
```

## Licence

Copyright © 2020–2026 the psvpp-ba authors, listed in [AUTHORS](AUTHORS).

* **Code** (everything outside `instances/`): GNU General Public License,
  version 3 or (at your option) any later version — see [LICENSE](LICENSE).
  Every source file carries the SPDX tags `SPDX-FileCopyrightText` and
  `SPDX-License-Identifier: GPL-3.0-or-later`.
* **Benchmark instances** (`instances/`): Creative Commons Attribution 4.0
  International (CC BY 4.0) — see [instances/LICENSE](instances/LICENSE). Their
  origin and how to credit them are in
  [instances/README.md](instances/README.md).
