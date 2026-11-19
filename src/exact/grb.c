// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file grb.c
 * @brief Gurobi helpers (see grb.h).
 */
#include "exact/grb.h"

#include "common/util.h"

#include <stdlib.h>

void grb_check(GRBenv *env, int error, const char *what, const char *file, int line)
{
    if (error == 0)
        return;
    die("Gurobi error %d in %s (%s:%d): %s", error, what, file, line,
        env ? GRBgeterrormsg(env) : "no environment available");
}

GRBenv *grb_open_env(bool log)
{
    GRBenv *env = nullptr;
    int err = GRBemptyenv(&env);
    if (err == 0)
        err = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, log ? 1 : 0);
    if (err == 0)
        err = GRBstartenv(env);
    if (err != 0)
        die("cannot start Gurobi (error %d): %s\n"
            "       check that a valid licence is installed (GRB_LICENSE_FILE or ~/gurobi.lic)",
            err, env ? GRBgeterrormsg(env) : "");
    return env;
}

GRBmodel *grb_new_model(GRBenv *env, const char *name, int threads, bool log)
{
    GRBmodel *m = nullptr;
    GRB_CHECK(env, GRBnewmodel(env, &m, name, 0, nullptr, nullptr, nullptr, nullptr, nullptr));
    GRBenv *menv = GRBgetenv(m);
    GRB_CHECK(env, GRBsetintattr(m, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE));
    /* the log is switched first: a silent model must not report its parameter changes */
    GRB_CHECK(env, GRBsetintparam(menv, GRB_INT_PAR_OUTPUTFLAG, log ? 1 : 0));
    GRB_CHECK(env, GRBsetintparam(menv, GRB_INT_PAR_THREADS, threads));
    return m;
}

void row_add(GrbRow *r, int var, double coef)
{
    if (r->len == r->cap) {
        r->cap = r->cap ? 2 * r->cap : 64;
        r->ind = xrealloc(r->ind, sizeof(int) * (size_t)r->cap);
        r->val = xrealloc(r->val, sizeof(double) * (size_t)r->cap);
    }
    r->ind[r->len] = var;
    r->val[r->len++] = coef;
}

void row_free(GrbRow *r)
{
    free(r->ind);
    free(r->val);
    *r = (GrbRow){0};
}

void grb_add_row(GRBmodel *model, const GrbRow *row, char sense, double rhs, const char *name)
{
    GRB_CHECK(GRBgetenv(model), GRBaddconstr(model, row->len, row->ind, row->val, sense, rhs, name));
}
