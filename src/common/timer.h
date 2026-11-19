// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file timer.h
 * @brief Wall-clock and CPU-time measurement.
 *
 * The heuristic reports the CPU time of the thread that executed the run, so
 * the figure stays meaningful when several runs are executed in parallel. The
 * exact methods report the CPU time of the whole process (all solver threads),
 * which is how the paper reports them, together with the wall-clock time.
 */
#ifndef PSVPP_COMMON_TIMER_H
#define PSVPP_COMMON_TIMER_H

/** @brief Seconds elapsed on a monotonic clock (origin unspecified). */
double wall_time(void);

/** @brief CPU seconds consumed by the calling thread. */
double thread_cpu_time(void);

/** @brief CPU seconds consumed by the whole process (all threads). */
double process_cpu_time(void);

#endif /* PSVPP_COMMON_TIMER_H */
