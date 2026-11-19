// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file util.h
 * @brief Small, dependency-free helpers shared by every module: fatal-error
 *        reporting, checked memory allocation and integer helpers.
 */
#ifndef PSVPP_COMMON_UTIL_H
#define PSVPP_COMMON_UTIL_H

#include <stddef.h>

/**
 * @brief Prints a formatted message to @c stderr (prefixed with "error: ")
 *        and terminates the program with @c EXIT_FAILURE.
 *
 * @param fmt printf-style format string followed by its arguments.
 * @return Never returns.
 */
[[noreturn, gnu::format(printf, 1, 2)]] void die(const char *fmt, ...);

/**
 * @brief @c malloc that never returns @c nullptr.
 *
 * @param size number of bytes (a request of 0 bytes allocates 1 byte so that
 *             the result is always a unique, freeable pointer).
 * @return Pointer to uninitialised memory. Aborts through die() on failure.
 */
[[nodiscard]] void *xmalloc(size_t size);

/**
 * @brief @c calloc that never returns @c nullptr.
 *
 * @param count number of elements.
 * @param size  size of each element in bytes.
 * @return Pointer to zero-initialised memory. Aborts through die() on failure
 *         (including multiplication overflow, which @c calloc detects).
 */
[[nodiscard]] void *xcalloc(size_t count, size_t size);

/**
 * @brief @c realloc that never returns @c nullptr.
 *
 * @param ptr  block to resize (may be @c nullptr).
 * @param size new size in bytes.
 * @return Pointer to the resized block. Aborts through die() on failure.
 */
[[nodiscard]] void *xrealloc(void *ptr, size_t size);

/**
 * @brief Duplicates a NUL-terminated string with xmalloc().
 *
 * @param s string to copy.
 * @return Newly allocated copy; release it with @c free.
 */
[[nodiscard]] char *xstrdup(const char *s);

/** @brief Allocates an uninitialised array of @p n elements of type @p T. */
#define XMALLOC_ARRAY(T, n) ((T *)xmalloc(sizeof(T) * (size_t)(n)))

/** @brief Allocates a zero-initialised array of @p n elements of type @p T. */
#define XCALLOC_ARRAY(T, n) ((T *)xcalloc((size_t)(n), sizeof(T)))

/** @brief Smaller of two integers. */
static inline int imin(int a, int b) { return a < b ? a : b; }

/** @brief Larger of two integers. */
static inline int imax(int a, int b) { return a > b ? a : b; }

/** @brief Positive part of an integer, i.e. max(x, 0). */
static inline int ipos(int x) { return x > 0 ? x : 0; }

/** @brief Absolute value of an integer. */
static inline int iabs(int x) { return x < 0 ? -x : x; }

/**
 * @brief Change of the positive part of @p base when @p delta is added to it:
 *        max(base + delta, 0) - max(base, 0).
 *
 * This is the elementary operation of every incremental penalty computation:
 * @p base is the current (possibly negative) slack violation of a constraint
 * and the result is how much the violation grows (or shrinks) after a move.
 */
static inline int excess_delta(int base, int delta) { return ipos(base + delta) - ipos(base); }

#endif /* PSVPP_COMMON_UTIL_H */
