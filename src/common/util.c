// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file util.c
 * @brief Implementation of the helpers declared in util.h.
 */
#include "common/util.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void die(const char *fmt, ...)
{
    fflush(stdout);
    fputs("error: ", stderr);
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

void *xmalloc(size_t size)
{
    void *p = malloc(size ? size : 1);
    if (p == nullptr)
        die("out of memory (requested %zu bytes)", size);
    return p;
}

void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count ? count : 1, size ? size : 1);
    if (p == nullptr)
        die("out of memory (requested %zu x %zu bytes)", count, size);
    return p;
}

void *xrealloc(void *ptr, size_t size)
{
    void *p = realloc(ptr, size ? size : 1);
    if (p == nullptr)
        die("out of memory (requested %zu bytes)", size);
    return p;
}

char *xstrdup(const char *s)
{
    const size_t n = strlen(s) + 1;
    char *copy = xmalloc(n);
    memcpy(copy, s, n);
    return copy;
}
