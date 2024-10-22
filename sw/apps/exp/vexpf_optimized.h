// Copyright 2024 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Luca Colagrande <colluca@iis.ee.ethz.ch>

#define BATCH_SIZE 8

#define ALLOCATE_BUFFER(type, size) \
    (type *)snrt_l1_alloc_cluster_local(size * sizeof(type), sizeof(type))

#define N_T_BUFFERS 2
#define N_KI_KD_BUFFERS 3

#include "vexpf_optimized_asm.h"

static inline void vexpf_optimized(float *a, float *b) {

    // Allocate input and output arrays of type double
    double *input = ALLOCATE_BUFFER(double, LEN);
    double *output = ALLOCATE_BUFFER(double, LEN);

    // Convert inputs to double
    for (int i = 0; i < LEN; i++) {
        input[i] = (double)a[i];
    }

    // Derived parameters
    unsigned int n_batches = LEN / BATCH_SIZE;
    unsigned int n_iterations = n_batches + 2;

    // Allocate buffers
    uint64_t *ki_buffers[N_KI_KD_BUFFERS];
    double *kd_buffers[N_KI_KD_BUFFERS];
    double *z_buffers[N_KI_KD_BUFFERS];
    uint64_t *t_buffers[N_T_BUFFERS];
    ki_buffers[0] = ALLOCATE_BUFFER(uint64_t, BATCH_SIZE);
    ki_buffers[1] = ALLOCATE_BUFFER(uint64_t, BATCH_SIZE);
    ki_buffers[2] = ALLOCATE_BUFFER(uint64_t, BATCH_SIZE);
    kd_buffers[0] = (double *)ki_buffers[0];
    kd_buffers[1] = (double *)ki_buffers[1];
    kd_buffers[2] = (double *)ki_buffers[2];
    z_buffers[0] = ALLOCATE_BUFFER(double, BATCH_SIZE);
    z_buffers[1] = ALLOCATE_BUFFER(double, BATCH_SIZE);
    z_buffers[2] = ALLOCATE_BUFFER(double, BATCH_SIZE);
    t_buffers[0] = ALLOCATE_BUFFER(uint64_t, BATCH_SIZE);
    t_buffers[1] = ALLOCATE_BUFFER(uint64_t, BATCH_SIZE);

    // Define buffer pointers for every phase (fp0, int and fp1)
    unsigned int fp0_k_idx = 0;
    unsigned int int_ki_idx = 0;
    unsigned int int_t_idx = 0;
    unsigned int fp1_kd_idx = 0;
    unsigned int fp1_t_idx = 0;
    double   *fp0_a_ptr = input;
    double   *fp0_k_ptr;
    double   *fp0_z_ptr;
    uint64_t *int_ki_ptr;
    uint64_t *int_t_ptr;
    double   *fp1_kd_ptr;
    uint64_t *fp1_t_ptr;
    double   *fp1_z_ptr;
    double   *fp1_b_ptr = output;

    // Exponential function constants
    uint32_t EXP2F_TABLE_BITS = 5;
    double N = 1 << EXP2F_TABLE_BITS;
    double InvLn2N = 0x1.71547652b82fep+0 * N;
    double SHIFT = 0x1.8p+52;
    double C[4] = {0x1.c6af84b912394p-5/N/N/N, 0x1.ebfce50fac4f3p-3/N/N, 0x1.62e42ff0c52d6p-1/N, 1.0};

    // Iterate over batches
    for (int iteration = 0; iteration < n_iterations; iteration++) {

        // FP0 phase
        if (iteration < n_iterations - 2) {

            // Index buffers
            fp0_k_ptr = kd_buffers[fp0_k_idx];
            fp0_z_ptr = z_buffers[fp0_k_idx];

            // Configure SSRs
            snrt_ssr_loop_1d(SNRT_SSR_DM_ALL, BATCH_SIZE, sizeof(double));
            snrt_ssr_read(SNRT_SSR_DM0, SNRT_SSR_1D, fp0_a_ptr);
            snrt_ssr_write(SNRT_SSR_DM1, SNRT_SSR_1D, fp0_k_ptr);
            snrt_ssr_write(SNRT_SSR_DM2, SNRT_SSR_1D, fp0_z_ptr);
            snrt_ssr_enable();

            // FP0 computation
            int unroll_factor = 4;
            for (int i = 0; i < BATCH_SIZE; i += unroll_factor) {
                asm volatile(
                    FP0_ASM_BODY
                    :
                    : [InvLn2N] "f" (InvLn2N), [SHIFT] "f" (SHIFT)
                    : "memory", "ft0", "ft1", "ft2", "fa3", "ft3", "ft4", "ft5"
                );
            }
            snrt_ssr_disable();

            // Increment input data pointer for next iteration
            fp0_a_ptr += BATCH_SIZE;

            // Increment buffer index for next iteration
            fp0_k_idx += 1;
            fp0_k_idx %= N_KI_KD_BUFFERS;
        }

        // INT phase
        if (iteration > 0 && iteration < n_iterations - 1) {

            // Index buffers
            int_ki_ptr = ki_buffers[int_ki_idx];
            int_t_ptr = t_buffers[int_t_idx];

            // INT computation
            int unroll_factor = 4;
            for (int i = 0; i < BATCH_SIZE; i += unroll_factor) {
                asm volatile(
                    INT_ASM_BODY
                    :
                    : [ki] "r" (int_ki_ptr + i), [T] "r" (T), [t] "r" (int_t_ptr + i)
                    : "memory", "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7",
                      "t0", "t1", "t2", "t3"
                );
            }

            // Increment buffer indices for next iteration
            int_ki_idx += 1;
            int_t_idx += 1;
            int_ki_idx %= N_KI_KD_BUFFERS;
            int_t_idx %= N_T_BUFFERS;
        }

        // FP1 phase
        if (iteration > 1) {

            // Index buffers
            fp1_kd_ptr = kd_buffers[fp1_kd_idx];
            fp1_t_ptr = t_buffers[fp1_t_idx];
            fp1_z_ptr = z_buffers[fp1_kd_idx];

            // FP1 computation
            int unroll_factor = 4;
            for (int i = 0; i < BATCH_SIZE; i += unroll_factor) {
                asm volatile(
                    FP1_ASM_BODY
                    :
                    : [kd] "r" (fp1_kd_ptr + i), [SHIFT] "f" (SHIFT),
                      [C0] "f" (C[0]), [C1] "f" (C[1]),
                      [C2] "f" (C[2]), [C3] "f" (C[3]),
                      [t] "r" (fp1_t_ptr + i), [output] "r" (fp1_b_ptr + i),
                      [z] "r" (fp1_z_ptr + i)
                    : "memory", "fa0", "fa1", "fa2", "fa3", "fa4", "fa5",
                      "fa6", "fa7", "ft3", "ft4", "ft5", "ft6", "ft7", "ft8",
                      "ft9", "ft10", "ft11", "fs0", "fs1", "fs2"
                );
            }

            // Increment input data pointer for next iteration
            fp1_b_ptr += BATCH_SIZE;

            // Increment buffer indices for next iteration
            fp1_kd_idx += 1;
            fp1_t_idx += 1;
            fp1_kd_idx %= N_KI_KD_BUFFERS;
            fp1_t_idx %= N_T_BUFFERS;
        }

        // Synchronize FP and integer threads
        snrt_fpu_fence();
    }

    // Convert outputs to float
    for (int i = 0; i < LEN; i++) {
        b[i] = (float)output[i];
    }
}
