// Copyright 2023 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Fabio Cappellini <fcappellini@student.ethz.ch>
// Hakim Filali <hfilali@student.ee.ethz.ch>
// Luca Colagrande <colluca@iis.ee.ethz.ch>
// Lannan Jiang <jiangl@student.ethz.ch>

#include "math.h"
#include "prng.h"
#include "snrt.h"

#ifndef N_SAMPLES
#define N_SAMPLES 1024
#endif

#ifndef FUNC_PTR
#define FUNC_PTR calculate_psum_naive
#endif

#ifndef N_CORES
#define N_CORES snrt_cluster_compute_core_num()
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 256
#endif

__thread double one = 1.0;

static inline uint32_t calculate_psum_naive(lcg_t *x_lcg, lcg_t *y_lcg,
                                            unsigned int n_samples) {
    uint32_t int_x = x_lcg->state;
    uint32_t int_y = y_lcg->state;
    double x, y;
    unsigned int result = 0;

    snrt_mcycle();
    for (unsigned int i = 0; i < n_samples; i++) {
        x = rand_int_to_unit_double(int_x);
        y = rand_int_to_unit_double(int_y);
        int_x = lcg_next(x_lcg);
        int_y = lcg_next(y_lcg);

        if ((x * x + y * y) < one) {
            result++;
        }
    }
    snrt_mcycle();

    return result;
}

static inline uint32_t calculate_psum_baseline(lcg_t *x_lcgs, lcg_t *y_lcgs,
                                               unsigned int n_samples) {
    if (snrt_cluster_core_idx() < N_CORES) {
        // Generate first 4 pseudo-random integer (X,Y) pairs
        unsigned int int_x0 = x_lcgs[0].state;
        unsigned int int_y0 = y_lcgs[0].state;
        unsigned int int_x1 = x_lcgs[1].state;
        unsigned int int_y1 = y_lcgs[1].state;
        unsigned int int_x2 = x_lcgs[2].state;
        unsigned int int_y2 = y_lcgs[2].state;
        unsigned int int_x3 = x_lcgs[3].state;
        unsigned int int_y3 = y_lcgs[3].state;
        int result = 0;
        int n_iter = n_samples / 4;

        snrt_mcycle();
        asm volatile(
            "1:"

            // Generate next 4 pseudo-random integer (X,Y) pairs
            "mul %[int_x0], %[int_x0], %[Ap_x] \n"
            "mul %[int_y0], %[int_y0], %[Ap_y] \n"
            "mul %[int_x1], %[int_x1], %[Ap_x] \n"
            "mul %[int_y1], %[int_y1], %[Ap_y] \n"
            "add %[int_x0], %[int_x0], %[Cp_x] \n"
            "mul %[int_x2], %[int_x2], %[Ap_x] \n"
            "add %[int_y0], %[int_y0], %[Cp_y] \n"
            "mul %[int_y2], %[int_y2], %[Ap_y] \n"
            "add %[int_x1], %[int_x1], %[Cp_x] \n"
            "mul %[int_x3], %[int_x3], %[Ap_x] \n"
            "add %[int_y1], %[int_y1], %[Cp_y] \n"
            "mul %[int_y3], %[int_y3], %[Ap_y] \n"
            "add %[int_x2], %[int_x2], %[Cp_x] \n"
            "add %[int_y2], %[int_y2], %[Cp_y] \n"
            "add %[int_x3], %[int_x3], %[Cp_x] \n"
            "add %[int_y3], %[int_y3], %[Cp_y] \n"

            // Convert integer PRNs to doubles
            "fcvt.d.wu ft0, %[int_x0] \n"
            "fcvt.d.wu fa0, %[int_y0] \n"
            "fcvt.d.wu ft1, %[int_x1] \n"
            "fcvt.d.wu fa1, %[int_y1] \n"
            "fcvt.d.wu ft2, %[int_x2] \n"
            "fcvt.d.wu fa2, %[int_y2] \n"
            "fcvt.d.wu ft3, %[int_x3] \n"
            "fcvt.d.wu fa3, %[int_y3] \n"

            // Normalize PRNs to [0, 1] range
            "fmul.d ft0, ft0, %[div] \n"
            "fmul.d ft1, ft1, %[div] \n"
            "fmul.d ft2, ft2, %[div] \n"
            "fmul.d ft3, ft3, %[div] \n"
            "fmul.d fa0, fa0, %[div] \n"
            "fmul.d fa1, fa1, %[div] \n"
            "fmul.d fa2, fa2, %[div] \n"
            "fmul.d fa3, fa3, %[div] \n"

            // X^2
            "fmul.d ft0, ft0, ft0 \n"
            "fmul.d ft1, ft1, ft1 \n"
            "fmul.d ft2, ft2, ft2 \n"
            "fmul.d ft3, ft3, ft3 \n"

            // X^2 + Y^2
            "fmadd.d ft0, fa0, fa0, ft0 \n"
            "fmadd.d ft1, fa1, fa1, ft1 \n"
            "fmadd.d ft2, fa2, fa2, ft2 \n"
            "fmadd.d ft3, fa3, fa3, ft3 \n"

            // (X^2 + Y^2) < 1
            "flt.d a0, ft0, %[one] \n"
            "flt.d a1, ft1, %[one] \n"
            "flt.d a2, ft2, %[one] \n"
            "flt.d a3, ft3, %[one] \n"

            // Count points in circle
            "add %[result], %[result], a0 \n"
            "add %[result], %[result], a1 \n"
            "add %[result], %[result], a2 \n"
            "add %[result], %[result], a3 \n"

            // Loop over batches
            "addi %[n_iter], %[n_iter], -1 \n"
            "bnez %[n_iter], 1b            \n"

            : [ result ] "+r"(result), [ n_iter ] "+r"(n_iter)
            : [ int_x0 ] "r"(int_x0), [ int_y0 ] "r"(int_y0),
              [ int_x1 ] "r"(int_x1), [ int_y1 ] "r"(int_y1),
              [ int_x2 ] "r"(int_x2), [ int_y2 ] "r"(int_y2),
              [ int_x3 ] "r"(int_x3), [ int_y3 ] "r"(int_y3),
              [ div ] "f"(max_uint_plus_1_inverse), [ one ] "f"(one),
              [ Ap_x ] "r"(x_lcgs->A), [ Cp_x ] "r"(x_lcgs->C),
              [ Ap_y ] "r"(y_lcgs->A), [ Cp_y ] "r"(y_lcgs->C)
            : "ft0", "ft1", "ft2", "ft3", "fa0", "fa1", "fa2", "fa3", "a0",
              "a1", "a2", "a3", "memory");
        snrt_fpu_fence();
        snrt_mcycle();

        return result;
    }

    return 0;
}

static inline uint32_t calculate_psum_optimized(lcg_t *x_lcgs, lcg_t *y_lcgs,
                                                unsigned int n_samples) {
    // Derived parameters
    uint32_t n_batches = n_samples / BATCH_SIZE;
    uint32_t n_iterations = n_batches + 2;

    // Allocate memory on TCDM
    unsigned int result = 0;
    double *rng_x_all[2];
    double *rng_y_all[2];
    double *rng_z_all[2];

    // Allocate (double) buffers for communication between integer and FP
    // threads in memory
    rng_x_all[0] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));
    rng_y_all[0] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));
    rng_z_all[0] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));
    rng_x_all[1] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));
    rng_y_all[1] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));
    rng_z_all[1] = (double *)snrt_l1_alloc_cluster_local(
        BATCH_SIZE * N_CORES * sizeof(double), sizeof(double));

    // Point each core to its own section of the buffers
    if (snrt_is_compute_core()) {
        unsigned int offset = snrt_cluster_core_idx() * BATCH_SIZE;
        rng_x_all[0] += offset;
        rng_y_all[0] += offset;
        rng_z_all[0] += offset;
        rng_x_all[1] += offset;
        rng_y_all[1] += offset;
        rng_z_all[1] += offset;
    }

    // Clear buffers. This is necessary since the PRNG generates 32-bit
    // integers, but the fcvt.w.wu.ssr instructions expects 64-bit integers.
    // Zero-padding is sufficient to encode a 64-bit unsigned integer from
    // a 32-bit integer. Thus, if the buffers are zero-initialized, storing the
    // 32-bit integers with a 64-bit stride is sufficient.
    if (snrt_is_dm_core()) {
        snrt_dma_memset(rng_x_all[0], 0, sizeof(double) * BATCH_SIZE * N_CORES);
        snrt_dma_memset(rng_y_all[0], 0, sizeof(double) * BATCH_SIZE * N_CORES);
        snrt_dma_memset(rng_z_all[0], 0, sizeof(double) * BATCH_SIZE * N_CORES);
        snrt_dma_memset(rng_x_all[1], 0, sizeof(double) * BATCH_SIZE * N_CORES);
        snrt_dma_memset(rng_y_all[1], 0, sizeof(double) * BATCH_SIZE * N_CORES);
        snrt_dma_memset(rng_z_all[1], 0, sizeof(double) * BATCH_SIZE * N_CORES);
    }
    snrt_cluster_hw_barrier();

    // Pointers to current set of buffers
    unsigned int fp_xyz_idx = 0;
    unsigned int int_xy_idx = 0;
    unsigned int int_z_idx = 0;
    double *fp_x_ptr = rng_x_all[fp_xyz_idx];
    double *fp_y_ptr = rng_y_all[fp_xyz_idx];
    double *fp_z_ptr = rng_z_all[fp_xyz_idx];
    double *int_x_ptr = rng_x_all[int_xy_idx];
    double *int_y_ptr = rng_y_all[int_xy_idx];
    double *int_z_ptr = rng_z_all[int_z_idx];

    // Initialize 4 pseudo-random integer (X,Y) pairs
    unsigned int int_x0 = x_lcgs[0].state;
    unsigned int int_y0 = y_lcgs[0].state;
    unsigned int int_x1 = x_lcgs[1].state;
    unsigned int int_y1 = y_lcgs[1].state;
    unsigned int int_x2 = x_lcgs[2].state;
    unsigned int int_y2 = y_lcgs[2].state;
    unsigned int int_x3 = x_lcgs[3].state;
    unsigned int int_y3 = y_lcgs[3].state;

    // Accumulators for partial sums
    int temp0 = 0;
    int temp1 = 0;
    int temp2 = 0;
    int temp3 = 0;
    int temp4 = 0;
    int temp5 = 0;
    int temp6 = 0;
    int temp7 = 0;

    if (snrt_cluster_core_idx() < N_CORES) {
        // Set up SSRs
        snrt_ssr_loop_1d(SNRT_SSR_DM_ALL, BATCH_SIZE,
                         snrt_cluster_compute_core_num() * sizeof(double));

        // Batch iterations
        for (int iteration = 0; iteration < n_iterations; iteration++) {
            snrt_mcycle();

            // Floating-point thread works on all but first and last iterations
            if (iteration > 0 && iteration < n_iterations - 1) {
                // Switch buffers for floating-point thread
                fp_xyz_idx ^= 1;
                fp_x_ptr = rng_x_all[fp_xyz_idx];
                fp_y_ptr = rng_y_all[fp_xyz_idx];
                fp_z_ptr = rng_z_all[fp_xyz_idx];

                // Point SSRs to current buffers
                snrt_ssr_read(SNRT_SSR_DM0, SNRT_SSR_1D, fp_x_ptr);
                snrt_ssr_read(SNRT_SSR_DM1, SNRT_SSR_1D, fp_y_ptr);
                snrt_ssr_write(SNRT_SSR_DM2, SNRT_SSR_1D, fp_z_ptr);

                // Fix register used by 1.0 constant
                register double reg_one asm("ft8") = one;

                // Enable SSRs
                snrt_ssr_enable();

                // Floating-point thread
                asm volatile(
                    // Unrolled by 2
                    "frep.o %[n_frep], 14, 0, 0 \n"

                    // Convert integer PRNs to doubles
                    // fcvt.d.wu.ssr fa1, ft0
                    // fcvt.d.wu.ssr fa3, ft0
                    // fcvt.d.wu.ssr fa2, ft1
                    // fcvt.d.wu.ssr fa4, ft1
                    ".word %[fcvt0] \n"
                    ".word %[fcvt1] \n"
                    ".word %[fcvt2] \n"
                    ".word %[fcvt3] \n"

                    // Normalize PRNs to [0, 1] range
                    "fmul.d fa1, fa1, %[div] \n"
                    "fmul.d fa3, fa3, %[div] \n"
                    "fmul.d fa2, fa2, %[div] \n"
                    "fmul.d fa4, fa4, %[div] \n"

                    // X^2
                    "fmul.d fa1, fa1, fa1 \n"
                    "fmul.d fa3, fa3, fa3 \n"

                    // X^2 + Y^2
                    "fmadd.d fa2, fa2, fa2, fa1 \n"
                    "fmadd.d fa4, fa4, fa4, fa3 \n"

                    // (X^2 + Y^2) < 1
                    // flt.d.ssr ft2, fa2, ft8
                    // flt.d.ssr ft2, fa4, ft8
                    ".word %[flt0] \n"
                    ".word %[flt1] \n"
                    :
                    : [ n_frep ] "r"(BATCH_SIZE / 2 - 1),
                      [ fcvt0 ] "i"(FCVT_D_WU_SSR(11, 0)),
                      [ fcvt1 ] "i"(FCVT_D_WU_SSR(13, 0)),
                      [ fcvt2 ] "i"(FCVT_D_WU_SSR(12, 1)),
                      [ fcvt3 ] "i"(FCVT_D_WU_SSR(14, 1)),
                      [ div ] "f"(max_uint_plus_1_inverse),
                      [ one ] "f"(reg_one), [ flt0 ] "i"(FLT_D_SSR(2, 12, 28)),
                      [ flt1 ] "i"(FLT_D_SSR(2, 14, 28))
                    : "ft0", "ft1", "ft2", "ft8", "fa1", "fa2", "fa3", "fa4",
                      "memory");
            }

            // Integer thread produces PRNs in all but last two iterations
            if (iteration < n_iterations - 2) {
                // Switch X, Y buffers for integer generation thread
                int_xy_idx ^= 1;
                int_x_ptr = rng_x_all[int_xy_idx];
                int_y_ptr = rng_y_all[int_xy_idx];

                // Fix register used by 1.0 constant
                uint32_t Ap_x = x_lcgs->A;
                uint32_t Ap_y = y_lcgs->A;
                uint32_t Cp_x = x_lcgs->C;
                uint32_t Cp_y = y_lcgs->C;

                // Unrolled by 4
                for (int j = 0; j < BATCH_SIZE; j += 4) {
                    asm volatile(
                        // Compute 4 integer PRN (X, Y) pairs for iteration i+2
                        // TODO(colluca): why can't the store be done after the
                        // computation?
                        "mul %[int_x0], %[int_x0], %[Ap_x] \n"
                        "mul %[int_y0], %[int_y0], %[Ap_y] \n"
                        "mul %[int_x1], %[int_x1], %[Ap_x] \n"
                        "mul %[int_y1], %[int_y1], %[Ap_y] \n"
                        "add %[int_x0], %[int_x0], %[Cp_x] \n"
                        "mul %[int_x2], %[int_x2], %[Ap_x] \n"
                        "add %[int_y0], %[int_y0], %[Cp_y] \n"
                        "mul %[int_y2], %[int_y2], %[Ap_y] \n"
                        "add %[int_x1], %[int_x1], %[Cp_x] \n"
                        "mul %[int_x3], %[int_x3], %[Ap_x] \n"
                        "add %[int_y1], %[int_y1], %[Cp_y] \n"
                        "mul %[int_y3], %[int_y3], %[Ap_y] \n"
                        "add %[int_x2], %[int_x2], %[Cp_x] \n"
                        "add %[int_y2], %[int_y2], %[Cp_y] \n"
                        "add %[int_x3], %[int_x3], %[Cp_x] \n"
                        "add %[int_y3], %[int_y3], %[Cp_y] \n"

                        // Store 4 integer PRN (X, Y) pairs for iteration i+1,
                        // zero-padded to 64-bit unsigned integers
                        "sw %[int_x0],  0(%[rng_x]) \n"
                        "sw %[int_y0],  0(%[rng_y]) \n"
                        "sw %[int_x1],  8(%[rng_x]) \n"
                        "sw %[int_y1],  8(%[rng_y]) \n"
                        "sw %[int_x2], 16(%[rng_x]) \n"
                        "sw %[int_y2], 16(%[rng_y]) \n"
                        "sw %[int_x3], 24(%[rng_x]) \n"
                        "sw %[int_y3], 24(%[rng_y]) \n"
                        : [ int_x0 ] "+r"(int_x0), [ int_y0 ] "+r"(int_y0),
                          [ int_x1 ] "+r"(int_x1), [ int_y1 ] "+r"(int_y1),
                          [ int_x2 ] "+r"(int_x2), [ int_y2 ] "+r"(int_y2),
                          [ int_x3 ] "+r"(int_x3), [ int_y3 ] "+r"(int_y3)
                        : [ rng_x ] "r"(&int_x_ptr[j]),
                          [ rng_y ] "r"(&int_y_ptr[j]), [ Ap_x ] "r"(Ap_x),
                          [ Cp_x ] "r"(Cp_x), [ Ap_y ] "r"(Ap_y),
                          [ Cp_y ] "r"(Cp_y)
                        : "memory");
                }
            }

            // Integer thread accumulates the comparison results in all
            // iterations after the second
            if (iteration > 1) {
                // Switch Z buffers for integer accumulation thread
                int_z_idx ^= 1;
                int_z_ptr = rng_z_all[int_z_idx];

                // Unrolled by 8
                for (int j = 0; j < BATCH_SIZE; j += 8) {
                    asm volatile(
                        "lw  a0,  0(%[rng_z])       \n"
                        "lw  a1,  8(%[rng_z])       \n"
                        "lw  a2, 16(%[rng_z])       \n"
                        "lw  a3, 24(%[rng_z])       \n"
                        "add %[temp0], %[temp0], a0 \n"
                        "lw  a4, 32(%[rng_z])       \n"
                        "add %[temp1], %[temp1], a1 \n"
                        "lw  a5, 40(%[rng_z])       \n"
                        "add %[temp2], %[temp2], a2 \n"
                        "lw  a6, 48(%[rng_z])       \n"
                        "add %[temp3], %[temp3], a3 \n"
                        "lw  a7, 56(%[rng_z])       \n"
                        "add %[temp4], %[temp4], a4 \n"
                        "add %[temp5], %[temp5], a5 \n"
                        "add %[temp6], %[temp6], a6 \n"
                        "add %[temp7], %[temp7], a7 \n"
                        : [ temp0 ] "+r"(temp0), [ temp1 ] "+r"(temp1),
                          [ temp2 ] "+r"(temp2), [ temp3 ] "+r"(temp3),
                          [ temp4 ] "+r"(temp4), [ temp5 ] "+r"(temp5),
                          [ temp6 ] "+r"(temp6), [ temp7 ] "+r"(temp7)
                        : [ rng_z ] "r"(&int_z_ptr[j])
                        : "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7",
                          "memory");
                }
            }

            // Synchronize FP and integer threads and disable SSRs
            snrt_fpu_fence();
            snrt_ssr_disable();
        }

        // Reduce partial sums
        result += temp0;
        result += temp1;
        result += temp2;
        result += temp3;
        result += temp4;
        result += temp5;
        result += temp6;
        result += temp7;
        return result;
    }

    return 0;
}

int main() {
    if (snrt_is_compute_core()) snrt_mcycle();

    // Initialize the PRNGs for parallel Monte Carlo.
    // Two independent X and Y sequences, divided into a subsequence per core.
    uint32_t n_seq = N_CORES;
    // In baseline implementation, each subsequence is further divided into 4
    // subsequences for unrolling.
    if (FUNC_PTR != calculate_psum_naive) n_seq *= 4;
    lcg_t *x_lcgs = (lcg_t *)snrt_l1_alloc_cluster_local(sizeof(lcg_t) * n_seq,
                                                         sizeof(lcg_t));
    lcg_t *y_lcgs = (lcg_t *)snrt_l1_alloc_cluster_local(sizeof(lcg_t) * n_seq,
                                                         sizeof(lcg_t));
    if (snrt_cluster_core_idx() == 0) {
        lcg_init_n_default(0, n_seq, x_lcgs);
        lcg_init_n_default(1, n_seq, y_lcgs);
    }
    snrt_cluster_hw_barrier();

    // Store partial sum array at first free address in TCDM
    uint32_t *reduction_array = (uint32_t *)snrt_l1_alloc_cluster_local(
        sizeof(uint32_t) * N_CORES, sizeof(uint32_t));

    // Calculate partial sums
    uint32_t n_samples_per_core = N_SAMPLES / N_CORES;
    uint32_t stride = (FUNC_PTR == calculate_psum_naive) ? 1 : 4;
    int result =
        FUNC_PTR(x_lcgs + snrt_cluster_core_idx() * stride,
                 y_lcgs + snrt_cluster_core_idx() * stride, n_samples_per_core);
    if (snrt_is_compute_core())
        reduction_array[snrt_cluster_core_idx()] = result;

    // Synchronize cores
    snrt_cluster_hw_barrier();

    // First core in cluster performs the final calculation
    if (snrt_cluster_core_idx() == 0) {
        // Reduce partial sums
        snrt_mcycle();
        uint32_t sum = 0;
        for (int i = 0; i < N_CORES; i++) {
            sum += reduction_array[i];
        }

        // Estimate pi
        double pi = (double)(4 * sum) / (double)N_SAMPLES;
        snrt_mcycle();

        // Check result
        double err = fabs(pi - M_PI);
        if (err > 0.1) return 1;
    }

    return 0;
}
