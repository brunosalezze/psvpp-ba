// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file timer.c
 * @brief POSIX implementation of the clocks declared in timer.h.
 */
/* clock_gettime() and its clocks are POSIX, not ISO C. */
#define _POSIX_C_SOURCE 200809L

#include "common/timer.h"

#include <time.h>

/* Seconds shown by POSIX clock @p id (0 if the clock is unavailable). */
static double read_clock(clockid_t id)
{
    struct timespec ts;
    if (clock_gettime(id, &ts) != 0)
        return 0.0;
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

double wall_time(void) { return read_clock(CLOCK_MONOTONIC); }

double thread_cpu_time(void) { return read_clock(CLOCK_THREAD_CPUTIME_ID); }

double process_cpu_time(void) { return read_clock(CLOCK_PROCESS_CPUTIME_ID); }
