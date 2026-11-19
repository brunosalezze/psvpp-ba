// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file rng.h
 * @brief Self-contained pseudo-random number generator (PCG32).
 *
 * Every heuristic run owns one ::Rng, so independent runs can be executed in
 * parallel and each run is reproducible from its seed alone, regardless of the
 * C library (the original code relied on the global, non-reentrant @c rand()).
 *
 * The generator is PCG-XSH-RR 64/32 by M. E. O'Neill (https://www.pcg-random.org).
 */
#ifndef PSVPP_COMMON_RNG_H
#define PSVPP_COMMON_RNG_H

#include <stdint.h>

/** @brief State of a PCG32 generator. Treat as opaque. */
typedef struct Rng {
    uint64_t state; /**< internal state */
    uint64_t inc;   /**< stream selector (always odd) */
} Rng;

/**
 * @brief Initialises a generator.
 *
 * @param rng    generator to initialise.
 * @param seed   any 64-bit value; equal seeds yield equal sequences.
 * @param stream stream identifier; generators with the same seed but different
 *               streams produce independent sequences.
 */
void rng_seed(Rng *rng, uint64_t seed, uint64_t stream);

/**
 * @brief Draws 32 uniformly distributed random bits.
 *
 * @param rng generator.
 * @return Next value of the sequence.
 */
uint32_t rng_u32(Rng *rng);

/**
 * @brief Draws an integer uniformly from the closed interval [@p lo, @p hi].
 *
 * Uses rejection sampling, so the result is unbiased.
 *
 * @param rng generator.
 * @param lo  smallest value that may be returned.
 * @param hi  largest value that may be returned; when @p hi < @p lo the
 *            function returns @p lo.
 * @return Value in [@p lo, @p hi].
 */
int rng_int(Rng *rng, int lo, int hi);

/**
 * @brief Draws a real number uniformly from [0, 1).
 *
 * @param rng generator.
 * @return Value with 32 random bits of resolution.
 */
double rng_unit(Rng *rng);

/**
 * @brief Shuffles an array in place (Fisher-Yates, every permutation equally likely).
 *
 * @param rng generator.
 * @param a   array to permute.
 * @param n   number of elements (nothing happens when @p n < 2).
 */
void rng_shuffle_int(Rng *rng, int *a, int n);

#endif /* PSVPP_COMMON_RNG_H */
