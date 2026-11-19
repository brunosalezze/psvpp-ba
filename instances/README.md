# Benchmark instances

The 30 instances used in the computational experiments of

> B. S. Vieira, G. M. Ribeiro, L. Bahiense, R. Cruz, A. B. Mendes, G. Laporte.
> Exact and heuristic algorithms for the fleet composition and periodic routing
> problem of offshore supply vessels with berth allocation decisions.
> European Journal of Operational Research 295(3):908–923, 2021.
> [doi:10.1016/j.ejor.2021.03.022](https://doi.org/10.1016/j.ejor.2021.03.022)

The naming scheme is `p<N>-<B>-<C>`: `N` installations, `B`
berths, `C` clusters.

## Origin

The instances are based on real operations of a Brazilian oil company in the
Campos Basin (paper, Section 5.1).

| instances | origin |
|-----------|--------|
| `p10-1-1`, `p15-1-1`, `p41-*`, `p79-*` | cases C10, C15, C41 and C79 of R. Cruz, A. B. Mendes, L. Bahiense, Y. Wu, *Integrating berth allocation decisions in a fleet composition and periodic routing problem of platform supply vessels*, European Journal of Operational Research 275(1):334–346, 2019, [doi:10.1016/j.ejor.2018.11.030](https://doi.org/10.1016/j.ejor.2018.11.030), including their berth variants and cluster arrangements |
| `p21-*`, `p30-*`, `p50-*`, `p60-*`, `p69-*` | created by Vieira et al. (2021) from C79 with 4 berths: offshore units removed at random down to the target number, coordinates redrawn inside the bounding box of C79's units |
| `p15-1-2` and the clustered variants of the created instances (`p21-1-2`, `p30-2-4`, `p50-3-7`, `p60-3-8`, `p69-4-8`) | cluster arrangements by Vieira et al. (2021), based on the units' locations |
| `p40-2-1`, `p66-4-1` | created by Vieira et al. (2021) in the same way, used only for parameter tuning |

The files are those used in the paper's experiments, unmodified.

## Licence

Copyright © 2019–2021 the authors of Cruz et al. (2019) and Vieira et al.
(2021). The instances are licensed under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/);
the full legal code is in [LICENSE](LICENSE). You may share and adapt them for
any purpose, provided you give appropriate credit, link to the licence and
indicate if you made changes. To give credit, cite Vieira et al. (2021) and,
for the C10, C15, C41 and C79 cases, Cruz et al. (2019).

The program code in this repository is licensed separately, under the GNU
General Public License v3.0 or later (see [../LICENSE](../LICENSE)).
