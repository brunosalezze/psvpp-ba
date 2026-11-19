// SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file rng.c
 * @brief PCG32 implementation (see rng.h).
 */
#include "common/rng.h"

void rng_seed(Rng *rng, uint64_t seed, uint64_t stream)
{
    rng->state = 0;
    rng->inc = (stream << 1u) | 1u;
    (void)rng_u32(rng);
    rng->state += seed;
    (void)rng_u32(rng);
}

uint32_t rng_u32(Rng *rng)
{
    const uint64_t old = rng->state;
    rng->state = old * 6364136223846793005ULL + rng->inc;
    const uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    const uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

int rng_int(Rng *rng, int lo, int hi)
{
    if (hi <= lo)
        return lo;
    const uint32_t bound = (uint32_t)((int64_t)hi - lo + 1);
    /* Reject the lowest (2^32 mod bound) values so that every residue is equally likely. */
    const uint32_t threshold = (uint32_t)(-bound) % bound;
    for (;;) {
        const uint32_t r = rng_u32(rng);
        if (r >= threshold)
            return lo + (int)(r % bound);
    }
}

double rng_unit(Rng *rng) { return (double)rng_u32(rng) * 0x1.0p-32; }

void rng_shuffle_int(Rng *rng, int *a, int n)
{
    for (int i = n - 1; i > 0; --i) {
        const int j = rng_int(rng, 0, i);
        const int t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
}
