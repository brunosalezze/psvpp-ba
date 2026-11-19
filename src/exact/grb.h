// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file grb.h
 * @brief Thin helpers around the Gurobi C API: error checking, model creation
 *        and a growable buffer for sparse rows.
 */
#ifndef PSVPP_EXACT_GRB_H
#define PSVPP_EXACT_GRB_H

#include <gurobi_c.h>

/**
 * @brief Aborts with Gurobi's error message when @p error is non-zero.
 *
 * @param env   environment that reports the message (may be @c nullptr).
 * @param error return code of a Gurobi call.
 * @param what  text describing the call (printed with the message).
 * @param file  source file of the call.
 * @param line  source line of the call.
 */
void grb_check(GRBenv *env, int error, const char *what, const char *file, int line);

/** @brief Calls a Gurobi function and aborts with a clear message on failure. */
#define GRB_CHECK(env, call) grb_check((env), (call), #call, __FILE__, __LINE__)

/**
 * @brief Creates the Gurobi environment (reads the licence).
 *
 * @param log true to let Gurobi print its log.
 * @return The environment; free it with @c GRBfreeenv.
 */
[[nodiscard]] GRBenv *grb_open_env(bool log);

/**
 * @brief Creates an empty minimisation model.
 *
 * @param env     environment.
 * @param name    model name.
 * @param threads Threads parameter (0 = all cores).
 * @param log     OutputFlag of the model.
 * @return The model; free it with @c GRBfreemodel.
 */
[[nodiscard]] GRBmodel *grb_new_model(GRBenv *env, const char *name, int threads, bool log);

/** @brief Sparse row (or column) under construction. */
typedef struct GrbRow {
    int *ind;    /**< variable indices */
    double *val; /**< coefficients */
    int len;     /**< entries in use */
    int cap;     /**< allocated entries */
} GrbRow;

/** @brief Empties a row, keeping its storage. */
static inline void row_clear(GrbRow *r) { r->len = 0; }

/** @brief Appends the term @p coef * x[@p var] to a row. */
void row_add(GrbRow *r, int var, double coef);

/** @brief Releases the storage of a row. */
void row_free(GrbRow *r);

/** @brief Adds @p row as a linear constraint `row sense rhs` to @p model. */
void grb_add_row(GRBmodel *model, const GrbRow *row, char sense, double rhs, const char *name);

#endif /* PSVPP_EXACT_GRB_H */
