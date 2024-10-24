// Copyright 2024 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Luca Colagrande <colluca@iis.ee.ethz.ch>

#include "math.h"
#include "snrt.h"

#ifndef LEN
#define LEN 12288
#endif

#ifndef FUNC_PTR
#define FUNC_PTR vexpf_optimized
#endif

__thread uint64_t T[64] = {
    0x3ff0000000000000, 0x3fefd9b0d3158574, 0x3fefb5586cf9890f, 0x3fef9301d0125b51,
    0x3fef72b83c7d517b, 0x3fef54873168b9aa, 0x3fef387a6e756238, 0x3fef1e9df51fdee1,
    0x3fef06fe0a31b715, 0x3feef1a7373aa9cb, 0x3feedea64c123422, 0x3feece086061892d,
    0x3feebfdad5362a27, 0x3feeb42b569d4f82, 0x3feeab07dd485429, 0x3feea47eb03a5585,
    0x3feea09e667f3bcd, 0x3fee9f75e8ec5f74, 0x3feea11473eb0187, 0x3feea589994cce13,
    0x3feeace5422aa0db, 0x3feeb737b0cdc5e5, 0x3feec49182a3f090, 0x3feed503b23e255d,
    0x3feee89f995ad3ad, 0x3feeff76f2fb5e47, 0x3fef199bdd85529c, 0x3fef3720dcef9069,
    0x3fef5818dcfba487, 0x3fef7c97337b9b5f, 0x3fefa4afa2a490da, 0x3fefd0765b6e4540,
};

double a[LEN], b_golden[LEN], b_actual[LEN];

#include "vexpf_naive.h"
#include "vexpf_baseline.h"
#include "vexpf_optimized.h"

int main() {
    uint32_t tstart, tend;

    // Initialize input array
    if (snrt_cluster_core_idx() == 0)
        for (int i = 0; i < LEN; i++)
            a[i] = (float)i / LEN;

    // Calculate exponential of input array using reference implementation
    if (snrt_cluster_core_idx() == 0) {
        tstart = snrt_mcycle();
        for (int i = 0; i < LEN; i++) {
            b_golden[i] = (double)expf((float)a[i]);
        }
        tend = snrt_mcycle();
    }

    // Synchronize cores
    snrt_cluster_hw_barrier();

    // Calculate exponential of input array using vectorized implementation
    tstart = snrt_mcycle();
    FUNC_PTR(a, b_actual);
    tend = snrt_mcycle();

    // Check if the results are correct
    if (snrt_cluster_core_idx() == 0) {
        uint32_t n_err = LEN;
        for (int i = 0; i < LEN; i++) {
            if ((float)b_golden[i] != (float)b_actual[i])
                printf("Error: b_golden[%d] = %f, b_actual[%d] = %f\n", i,
                    (float)b_golden[i], i, (float)b_actual[i]);
            else
                n_err--;
        }
        return n_err;
    } else
        return 0;
}
