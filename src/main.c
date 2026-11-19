// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file main.c
 * @brief Command-line interface of psvpp-ba.
 *
 * Usage: `psvpp-ba <command> [options] <instance>`; run `psvpp-ba --help`.
 *
 * Commands:
 *  - `alns`     multi-start ALNS heuristic (paper, Section 4);
 *  - `routes`   enumeration of the feasible voyages used by the exact models;
 *  - `bc`       branch-and-cut on the reduced model (Sections 2.2 and 3, needs Gurobi);
 *  - `complete` complete model (Section 2.1, needs Gurobi);
 *  - `info`     instance summary.
 */
#include "common/instance.h"
#include "common/timer.h"
#include "common/util.h"
#include "exact/exact.h"
#include "heuristic/alns.h"

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PSVPP_VERSION
#define PSVPP_VERSION "1.0.0"
#endif

/* ------------------------------------------------------------------------- */
/*  Help                                                                     */
/* ------------------------------------------------------------------------- */

/* Prints the help text, with the default value of every option, to @p out. */
static void print_usage(FILE *out)
{
    const AlnsParams d = alns_default_params();
    const ExactParams e = exact_default_params();
    fprintf(out,
            "psvpp-ba %s - fleet composition and periodic routing of offshore supply vessels\n"
            "           with berth allocation decisions (PSVPP-BA)\n\n"
            "usage: psvpp-ba <command> [options] <instance-file>\n\n"
            "commands:\n"
            "  alns       multi-start adaptive large neighbourhood search (heuristic)\n"
            "  routes     enumerate the feasible voyages used by the exact models\n"
            "  bc         branch-and-cut on the reduced model (requires Gurobi)\n"
            "  complete   solve the complete model (requires Gurobi)\n"
            "  info       print a summary of the instance\n\n"
            "general options:\n"
            "  -o, --output FILE        write the best solution to FILE\n"
            "  -q, --quiet              print only the result lines\n"
            "  -v, --verbose            alns: print the progress of the starts and the operator statistics\n"
            "      --max-route-size N   installations per voyage, at most (%d)\n"
            "  -h, --help               show this help\n"
            "  -V, --version            show the version\n\n"
            "alns options (defaults are the tuned values of the paper):\n"
            "  -s, --seed N             seed of the first run (1); run r uses seed N + r - 1\n"
            "  -r, --runs N             independent runs (1)\n"
            "  -j, --threads N          runs executed in parallel (1; needs OpenMP)\n"
            "      --csv                print one CSV line per run\n"
            "      --starts N           alpha: number of starts (%d)\n"
            "      --start-iterations N beta: iterations per start (%d)\n"
            "      --iterations N       gamma: total iterations (%d)\n"
            "      --ls-interval N      delta: iterations between local searches (%d)\n"
            "      --ls-max-rounds N    improvement steps per local search, at most (%d)\n"
            "      --segment N          rho: iterations per adaptive segment (%d)\n"
            "      --t0 X               initial temperature relative to the cost (%g)\n"
            "      --tf X               final temperature relative to the cost (%g)\n"
            "      --removal-min X      lower fraction of removed installations (%g)\n"
            "      --removal-max X      upper fraction of removed installations (%g)\n"
            "      --determinism X      p: determinism of the guided removals (%g)\n"
            "      --reaction X         tau_A: reaction factor of the adaptive layer (%g)\n"
            "      --scores A,B,C       sigma_1..3: scores of best/improving/accepted (%g,%g,%g)\n"
            "      --uniform            choose operators uniformly (disable the adaptive layer)\n"
            "      --penalty-increase X factor for too infeasible segments (%g)\n"
            "      --penalty-decrease X factor for too feasible segments (%g)\n\n"
            "exact options (routes, bc, complete):\n"
            "      --time-limit S       time limit in seconds (%g)\n"
            "      --gap PCT            stop when the gap is at most PCT%% (%g, Gurobi's default;\n"
            "                           the original code stopped the branch-and-cut at 0.1)\n"
            "      --solver-threads N   Gurobi threads of the main model, 0 = all cores (%d)\n"
            "      --subproblem-threads N\n"
            "                           bc: Gurobi threads of the fallback models of the separation,\n"
            "                           0 = all cores, at most --solver-threads (%d)\n"
            "      --min-usage PCT      bc: drop voyages loading less than PCT%% of the vessel on\n"
            "                           non-clustered instances; 0 keeps the method exact (%g)\n"
            "      --root-only          bc: only solve the LP relaxation and report its bound\n"
            "      --solver-log         show the Gurobi log\n",
            PSVPP_VERSION, e.max_route_size, d.starts, d.start_iterations, d.iterations, d.ls_interval, d.ls_max_rounds,
            d.segment, d.t0, d.tf, d.removal_min, d.removal_max, d.determinism, d.reaction, d.score_best,
            d.score_better, d.score_accept, d.penalty_increase, d.penalty_decrease, e.time_limit, 100 * e.gap,
            e.threads, e.subproblem_threads, e.min_usage);
}

/* ------------------------------------------------------------------------- */
/*  Option parsing                                                           */
/* ------------------------------------------------------------------------- */

typedef enum Command { CMD_ALNS, CMD_ROUTES, CMD_BC, CMD_COMPLETE, CMD_INFO } Command;

typedef struct Options {
    Command cmd;
    const char *instance;
    const char *output;
    bool quiet;
    bool verbose;
    bool csv;
    uint64_t seed;
    int runs;
    int threads;
    AlnsParams alns;
    ExactParams exact;
} Options;

enum {
    OPT_CSV = 256,
    OPT_LS_ROUNDS,
    OPT_STARTS,
    OPT_START_ITER,
    OPT_ITER,
    OPT_LS,
    OPT_SEGMENT,
    OPT_T0,
    OPT_TF,
    OPT_RMIN,
    OPT_RMAX,
    OPT_DET,
    OPT_REACTION,
    OPT_SCORES,
    OPT_UNIFORM,
    OPT_PINC,
    OPT_PDEC,
    OPT_TIME,
    OPT_SOLVER_THREADS,
    OPT_SUBPROBLEM_THREADS,
    OPT_ROUTE_SIZE,
    OPT_MIN_USAGE,
    OPT_GAP,
    OPT_ROOT_ONLY,
    OPT_SOLVER_LOG,
};

/* Parses @p arg as an integer in [min, max]; aborts naming option @p opt if it is not one. */
static long parse_long(const char *arg, const char *opt, long min, long max)
{
    char *end;
    errno = 0;
    const long v = strtol(arg, &end, 10);
    if (*arg == '\0' || *end != '\0' || errno != 0 || v < min || v > max)
        die("invalid value '%s' for %s (expected an integer in [%ld, %ld])", arg, opt, min, max);
    return v;
}

/* Parses @p arg as a real number in [min, max]; aborts naming option @p opt if it is not one. */
static double parse_double(const char *arg, const char *opt, double min, double max)
{
    char *end;
    const double v = strtod(arg, &end);
    if (*arg == '\0' || *end != '\0' || !(v >= min && v <= max))
        die("invalid value '%s' for %s (expected a number in [%g, %g])", arg, opt, min, max);
    return v;
}

/* Parses the command line into an ::Options value; prints help/version and exits when asked. */
static Options parse_options(int argc, char **argv)
{
    if (argc >= 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        print_usage(stdout);
        exit(EXIT_SUCCESS);
    }
    if (argc >= 2 && (!strcmp(argv[1], "-V") || !strcmp(argv[1], "--version"))) {
        printf("psvpp-ba %s\n", PSVPP_VERSION);
        exit(EXIT_SUCCESS);
    }
    if (argc < 2) {
        print_usage(stderr);
        exit(EXIT_FAILURE);
    }
    Options o = {.seed = 1, .runs = 1, .threads = 1, .alns = alns_default_params(), .exact = exact_default_params()};
    const char *cmd = argv[1];
    if (!strcmp(cmd, "alns"))
        o.cmd = CMD_ALNS;
    else if (!strcmp(cmd, "routes"))
        o.cmd = CMD_ROUTES;
    else if (!strcmp(cmd, "bc"))
        o.cmd = CMD_BC;
    else if (!strcmp(cmd, "complete"))
        o.cmd = CMD_COMPLETE;
    else if (!strcmp(cmd, "info"))
        o.cmd = CMD_INFO;
    else
        die("unknown command '%s' (run 'psvpp-ba --help')", cmd);

    static const struct option longopts[] = {
        {"output", required_argument, nullptr, 'o'},
        {"quiet", no_argument, nullptr, 'q'},
        {"verbose", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        {"seed", required_argument, nullptr, 's'},
        {"runs", required_argument, nullptr, 'r'},
        {"threads", required_argument, nullptr, 'j'},
        {"csv", no_argument, nullptr, OPT_CSV},
        {"starts", required_argument, nullptr, OPT_STARTS},
        {"start-iterations", required_argument, nullptr, OPT_START_ITER},
        {"iterations", required_argument, nullptr, OPT_ITER},
        {"ls-interval", required_argument, nullptr, OPT_LS},
        {"ls-max-rounds", required_argument, nullptr, OPT_LS_ROUNDS},
        {"segment", required_argument, nullptr, OPT_SEGMENT},
        {"t0", required_argument, nullptr, OPT_T0},
        {"tf", required_argument, nullptr, OPT_TF},
        {"removal-min", required_argument, nullptr, OPT_RMIN},
        {"removal-max", required_argument, nullptr, OPT_RMAX},
        {"determinism", required_argument, nullptr, OPT_DET},
        {"reaction", required_argument, nullptr, OPT_REACTION},
        {"scores", required_argument, nullptr, OPT_SCORES},
        {"uniform", no_argument, nullptr, OPT_UNIFORM},
        {"penalty-increase", required_argument, nullptr, OPT_PINC},
        {"penalty-decrease", required_argument, nullptr, OPT_PDEC},
        {"time-limit", required_argument, nullptr, OPT_TIME},
        {"solver-threads", required_argument, nullptr, OPT_SOLVER_THREADS},
        {"subproblem-threads", required_argument, nullptr, OPT_SUBPROBLEM_THREADS},
        {"max-route-size", required_argument, nullptr, OPT_ROUTE_SIZE},
        {"min-usage", required_argument, nullptr, OPT_MIN_USAGE},
        {"gap", required_argument, nullptr, OPT_GAP},
        {"root-only", no_argument, nullptr, OPT_ROOT_ONLY},
        {"solver-log", no_argument, nullptr, OPT_SOLVER_LOG},
        {nullptr, 0, nullptr, 0},
    };
    optind = 2;
    int c;
    while ((c = getopt_long(argc, argv, "o:qvhs:r:j:", longopts, nullptr)) != -1) {
        AlnsParams *a = &o.alns;
        switch (c) {
        case 'o': o.output = optarg; break;
        case 'q': o.quiet = true; break;
        case 'v': o.verbose = true; break;
        case 'h': print_usage(stdout); exit(EXIT_SUCCESS);
        case 's': o.seed = (uint64_t)parse_long(optarg, "--seed", 0, INT32_MAX); break;
        case 'r': o.runs = (int)parse_long(optarg, "--runs", 1, 1000000); break;
        case 'j': o.threads = (int)parse_long(optarg, "--threads", 1, 4096); break;
        case OPT_CSV: o.csv = true; break;
        case OPT_STARTS: a->starts = (int)parse_long(optarg, "--starts", 0, INT32_MAX); break;
        case OPT_START_ITER: a->start_iterations = (int)parse_long(optarg, "--start-iterations", 1, INT32_MAX); break;
        case OPT_ITER: a->iterations = (int)parse_long(optarg, "--iterations", 1, INT32_MAX); break;
        case OPT_LS: a->ls_interval = (int)parse_long(optarg, "--ls-interval", 1, INT32_MAX); break;
        case OPT_LS_ROUNDS: a->ls_max_rounds = (int)parse_long(optarg, "--ls-max-rounds", 0, INT32_MAX); break;
        case OPT_SEGMENT: a->segment = (int)parse_long(optarg, "--segment", 1, INT32_MAX); break;
        case OPT_T0: a->t0 = parse_double(optarg, "--t0", 1e-12, 1e6); break;
        case OPT_TF: a->tf = parse_double(optarg, "--tf", 1e-12, 1e6); break;
        case OPT_RMIN: a->removal_min = parse_double(optarg, "--removal-min", 0, 1); break;
        case OPT_RMAX: a->removal_max = parse_double(optarg, "--removal-max", 0, 1); break;
        case OPT_DET: a->determinism = parse_double(optarg, "--determinism", 0, 1e3); break;
        case OPT_REACTION: a->reaction = parse_double(optarg, "--reaction", 0, 1); break;
        case OPT_SCORES:
            if (sscanf(optarg, "%lf,%lf,%lf", &a->score_best, &a->score_better, &a->score_accept) != 3 ||
                !(a->score_best >= 0 && a->score_better >= 0 && a->score_accept >= 0))
                die("invalid value '%s' for --scores (expected three non-negative numbers: A,B,C)", optarg);
            break;
        case OPT_UNIFORM: a->adaptive = false; break;
        case OPT_PINC: a->penalty_increase = parse_double(optarg, "--penalty-increase", 1, 100); break;
        case OPT_PDEC: a->penalty_decrease = parse_double(optarg, "--penalty-decrease", 0, 1); break;
        case OPT_TIME: o.exact.time_limit = parse_double(optarg, "--time-limit", 0, 1e9); break;
        case OPT_SOLVER_THREADS: o.exact.threads = (int)parse_long(optarg, "--solver-threads", 0, 4096); break;
        case OPT_SUBPROBLEM_THREADS:
            o.exact.subproblem_threads = (int)parse_long(optarg, "--subproblem-threads", 0, 4096);
            break;
        case OPT_ROUTE_SIZE:
            o.exact.max_route_size = (int)parse_long(optarg, "--max-route-size", 1, MAX_ROUTE_SIZE_LIMIT);
            a->max_route_size = o.exact.max_route_size;
            break;
        case OPT_MIN_USAGE: o.exact.min_usage = parse_double(optarg, "--min-usage", 0, 100); break;
        case OPT_GAP: o.exact.gap = parse_double(optarg, "--gap", 0, 100) / 100.0; break;
        case OPT_SOLVER_LOG: o.exact.solver_log = true; break;
        case OPT_ROOT_ONLY: o.exact.root_only = true; break;
        default: die("run 'psvpp-ba --help' for the list of options");
        }
    }
    if (optind != argc - 1)
        die("expected exactly one instance file after the options (run 'psvpp-ba --help')");
    o.instance = argv[optind];
    if (o.csv)
        o.verbose = false; /* the progress messages would break the CSV output */
    if (o.alns.removal_min > o.alns.removal_max)
        die("--removal-min must not exceed --removal-max");
    return o;
}

/* ------------------------------------------------------------------------- */
/*  alns command                                                             */
/* ------------------------------------------------------------------------- */

typedef struct RunResult {
    uint64_t seed;
    double cost;      /* penalised objective, reported units */
    double objective; /* true objective (fixed + sailing), reported units */
    bool feasible;
    double cpu; /* CPU seconds of the thread that executed the run */
} RunResult;

/* Prints the result line of run @p r (or its CSV record): objective, feasibility, fleet and CPU time. */
static void print_run(const Options *o, const Instance *inst, int r, const RunResult *res, const Solution *s)
{
    int *fleet = XMALLOC_ARRAY(int, inst->num_types);
    sol_fleet(s, fleet);
    if (o->csv) {
        printf("%s,%llu,%.6f,%.6f,%d,%.3f", inst->name, (unsigned long long)res->seed, res->objective, res->cost,
               res->feasible, res->cpu);
        for (int k = 0; k < inst->num_types; k++)
            printf(",%d", fleet[k]);
        printf("\n");
    } else {
        printf("run %d (seed %llu): objective %.6f%s, fleet [", r + 1, (unsigned long long)res->seed, res->objective,
               res->feasible ? "" : " INFEASIBLE");
        for (int k = 0; k < inst->num_types; k++)
            printf("%s%d", k ? " " : "", fleet[k]);
        printf("], cpu %.2f s\n", res->cpu);
    }
    fflush(stdout);
    free(fleet);
}

/* Prints the statistics of a run: iterations, chosen start, and uses and final weight of each operator. */
static void print_report(const AlnsReport *rep)
{
    printf("  %ld iterations (%ld with a feasible candidate); final ALNS started from start %d (%.6f)\n",
           rep->iterations, rep->feasible_iterations, rep->best_start + 1, rep->start_cost / COST_SCALE);
    printf("  operator          uses   final weight\n");
    for (int i = 0; i < NUM_REMOVAL_OPS; i++)
        printf("  %-18s %8ld   %.3f\n", removal_op_name((RemovalOp)i), rep->removal_uses[i], rep->removal_weight[i]);
    for (int i = 0; i < NUM_INSERTION_OPS; i++)
        printf("  %-18s %8ld   %.3f\n", insertion_op_name((InsertionOp)i), rep->insertion_uses[i],
               rep->insertion_weight[i]);
}

/* `alns` command: performs the requested runs (in parallel with OpenMP), prints each result, the
 * summary statistics and the best solution, and writes it to the output file. Returns the exit status. */
static int run_alns(const Options *o, const Instance *inst)
{
    RunResult *res = XCALLOC_ARRAY(RunResult, o->runs);
    Solution *best = sol_create(inst);
    bool have_best = false;
    double best_objective = INFINITY;
    const double wall0 = wall_time();
    if (!o->quiet && !o->csv)
        printf("ALNS: %d run(s), alpha=%d beta=%d gamma=%d delta=%d rho=%d, %s operator selection\n", o->runs,
               o->alns.starts, o->alns.start_iterations, o->alns.iterations, o->alns.ls_interval, o->alns.segment,
               o->alns.adaptive ? "adaptive" : "uniform");
    if (o->csv) {
        printf("instance,seed,objective,penalised,feasible,cpu_s");
        for (int k = 0; k < inst->num_types; k++)
            printf(",fleet_type%d", inst->type[k].id);
        printf("\n");
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(o->threads)
#endif
    for (int r = 0; r < o->runs; r++) {
        Search S;
        search_init(&S, inst, &o->alns, o->seed + (uint64_t)r);
        Solution *s = sol_create(inst);
        const double c0 = thread_cpu_time();
        AlnsReport report;
        /* progress messages only make sense when runs do not interleave */
        FILE *log = o->verbose && (o->runs == 1 || o->threads == 1) ? stdout : nullptr;
        alns_multistart(&S, s, &report, log);
        res[r] = (RunResult){.seed = o->seed + (uint64_t)r,
                             .cost = s->cost / COST_SCALE,
                             .objective = sol_objective(s),
                             .feasible = !s->infeasible,
                             .cpu = thread_cpu_time() - c0};
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            print_run(o, inst, r, &res[r], s);
            if (o->verbose)
                print_report(&report);
            const bool improves = s->infeasible != best->infeasible ? !s->infeasible
                                  : s->infeasible                   ? s->cost < best->cost
                                                                    : res[r].objective < best_objective;
            if (!have_best || improves) {
                sol_copy(best, s);
                best_objective = res[r].objective;
                have_best = true;
            }
        }
        sol_free(s);
        search_free(&S);
    }

    if (o->runs > 1 && !o->csv) {
        double sum = 0, sq = 0, cpu = 0, bestv = INFINITY;
        int feasible = 0;
        for (int r = 0; r < o->runs; r++) {
            sum += res[r].objective;
            sq += res[r].objective * res[r].objective;
            cpu += res[r].cpu;
            feasible += res[r].feasible;
            if (res[r].feasible && res[r].objective < bestv)
                bestv = res[r].objective;
        }
        const double avg = sum / o->runs;
        const double var = o->runs > 1 ? fmax(0.0, (sq - o->runs * avg * avg) / (o->runs - 1)) : 0.0;
        printf("summary: best %.6f, average %.6f, deviation %.2f%%, feasible %d/%d, average cpu %.2f s, wall %.2f s\n",
               bestv, avg, avg > 0 ? 100.0 * sqrt(var) / avg : 0.0, feasible, o->runs, cpu / o->runs,
               wall_time() - wall0);
    }
    if (!o->quiet && !o->csv) {
        printf("\nbest solution:\n");
        sol_print(best, stdout);
    }
    if (o->output) {
        FILE *f = fopen(o->output, "w");
        if (f == nullptr)
            die("cannot write '%s'", o->output);
        sol_print(best, f);
        fclose(f);
    }
    sol_free(best);
    free(res);
    return EXIT_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*  main                                                                     */
/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    const Options o = parse_options(argc, argv);
#ifndef _OPENMP
    if (o.threads > 1)
        fprintf(stderr, "warning: built without OpenMP; runs are executed sequentially\n");
#endif
    Instance *inst = inst_load(o.instance);
    int status = EXIT_SUCCESS;
    switch (o.cmd) {
    case CMD_INFO: inst_print_summary(inst); break;
    case CMD_ALNS:
        if (!o.quiet && !o.csv)
            inst_print_summary(inst);
        status = run_alns(&o, inst);
        break;
    case CMD_ROUTES: status = exact_routes_command(inst, &o.exact, o.quiet); break;
    case CMD_BC: status = exact_bc_command(inst, &o.exact, o.quiet, o.output); break;
    case CMD_COMPLETE: status = exact_complete_command(inst, &o.exact, o.quiet, o.output); break;
    }
    inst_free(inst);
    return status;
}
