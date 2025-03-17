// Copyright 2023 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Author: Tim Fischer <fischeti@iis.ee.ethz.ch>
//         Luca Bertaccini <lbertaccini@iis.ee.ethz.ch>
//         Luca Colagrande <colluca@iis.ee.ethz.ch>
//         Viviane Potocnik <vivianep@iis.ee.ethz.ch>

#include <stdint.h>

#include "snrt.h"

#pragma once

typedef float v2f32 __attribute__((vector_size(8)));
typedef __fp16 v4f16 __attribute__((vector_size(8)));
typedef char v8f8 __attribute__((vector_size(8)));


#define ALIGN_NEXT_FROM_BASE(addr, base, size) \
    (((((addr) - (base)) + (size) - 1) / (size)) * (size) + (base))

#define BANK_WIDTH 8
#define HYPERBANK_WIDTH (SNRT_TCDM_BANK_PER_HYPERBANK_NUM * BANK_WIDTH)
#define ALIGN_UP_TCDM(addr) ALIGN_NEXT_FROM_BASE(addr, SNRT_TCDM_START_ADDR, HYPERBANK_WIDTH)

#include "gemm_fp16.h"
#include "gemm_fp32.h"
#include "gemm_fp64.h"
#include "gemm_fp8.h"

// define the gemm_fp function pointer
typedef void (*gemm_fp_t)(uint32_t m, uint32_t n, uint32_t k, void* a,
                          uint32_t lda, uint32_t transa, void* b,
                          uint32_t transb, uint32_t ldb, void* c, uint32_t ldc,
                          uint32_t beta, uint32_t setup_ssr);

typedef struct {
    double alpha;
    uint32_t prec;
    uint32_t setup_ssr;
    uint32_t parallelize_m;
    uint32_t parallelize_k;
    uint32_t m_tiles;
    uint32_t n_tiles;
    uint32_t k_tiles;
    uint32_t load_a;
    uint32_t load_b;
    uint32_t load_c;
    uint32_t transa;
    uint32_t transb;
    uint32_t M;
    uint32_t N;
    uint32_t K;
    void* a;
    void* b;
    uint32_t beta;
    void* c;
    void* gemm_fp;
} gemm_args_t;

// BLAS compliant single-cluster single-tile GEMM kernel, with some additional
// arguments at the beginning to specify Snitch implementation details. Matrix
// sizes and pointers are for the whole cluster computation. Within a cluster
// the computation is parallelized by assigning distinct output rows to
// distinct cores.
// TODO: beta (and alpha) should be of floating-point type (same precision as
// operands)
static inline void sc_st_gemm(gemm_args_t* gemm_args, void* a, void* b,
                              uint32_t beta, void* c) {
    gemm_fp_t impl = (gemm_fp_t)gemm_args->gemm_fp;
    precision_t prec = gemm_args->prec;
    uint32_t setup_ssr = gemm_args->setup_ssr;
    uint32_t transa = gemm_args->transa;
    uint32_t transb = gemm_args->transb;

    uint32_t m = gemm_args->M / gemm_args->m_tiles;
    uint32_t n = gemm_args->N / gemm_args->n_tiles;
    uint32_t k = gemm_args->K;

    uint32_t lda = k;
    uint32_t ldb;
    if (transb) {
        ldb = k;
    } else {
        ldb = n;
    }
    uint32_t ldc = n;
    int elements_per_line = 8;

    double alpha = gemm_args->alpha;

    if (snrt_is_compute_core()) {
        const uint32_t compute_num = snrt_cluster_compute_core_num();
        const uint32_t compute_id = snrt_cluster_core_idx();

        // Compute cores work not on contiguo usblocks but on strided rows
        uint32_t lda_strided = compute_num * lda;
        uint32_t ldc_strided = compute_num * ldc;

        //uint32_t lda_strided = compute_num * lda * 2;

        // Compute cores access A and C at offsets of one row from each other
        //uint32_t offsetA = compute_id * lda * prec;
        //uint32_t offsetC = compute_id * ldc * prec;

        // Compute cores access A and C at offsets of one row from each other Split Memory
        uint32_t offsetA = (compute_id * lda)/elements_per_line*HYPERBANK_WIDTH;
        uint32_t offsetC = (compute_id * ldc)/elements_per_line*HYPERBANK_WIDTH;

        //uint32_t offsetC;
        //if((compute_id*ldc)%16==0){
        //    offsetC = (compute_id * ldc)/16 * HYPERBANK_WIDTH;
        //}else{  //only case N=8
        //    offsetC = (compute_id * ldc)/16 * HYPERBANK_WIDTH+8*8 ;
        //}


        //uint32_t offsetA = compute_id * (int)(lda/16) * prec * 32;
        //uint32_t offsetA = (compute_id * lda * prec)/(HYPERBANK_WIDTH/2)*HYPERBANK_WIDTH+(compute_id * lda * prec)%(HYPERBANK_WIDTH/2);
        

        // Compute fraction of C rows every core computes
        uint32_t frac_m = m / compute_num;
        uint32_t rem_m = m % compute_num;
        if (snrt_cluster_core_idx() < rem_m) frac_m++;

        if (frac_m > 0)
            impl(frac_m, n, k, a + offsetA, lda_strided, transa, b, ldb, transb,
                 c + offsetC, ldc_strided, (float)beta, setup_ssr);
    }
}

// Multiple-cluster multiple-tile GEMM implementation.
// If parallelize_m, assigns a distinct subset of M-tiles to distinct clusters.
// If parallelize_k, then K-tiles are distributed to distinct clusters; a
// binary reduction tree is implemented to accumulate these tiles together.
// Note: in the current implementation, parallelize_m and parallelize_k
// should be mutually-exclusive. The load_* options allow to bypass the DMA
// transfers and operate directly on the a, b and c inputs.
// m_tiles: number of tiles in M dimension
// k_tiles: number of tiles in K dimension
// n_tiles: number of tiles in N dimension
static inline int gemm(const gemm_args_t* args) {

    // gemm_args_t* local_args = snrt_l1_next();
    gemm_args_t local_args;

    // // Copy the arguments to local memory
    // if (snrt_is_dm_core()) {
    //     snrt_dma_start_1d(local_args, args, sizeof(gemm_args_t));
    //     snrt_dma_wait_all();
    // }
    // snrt_cluster_hw_barrier();
    local_args = *args;

    uint32_t m = local_args.M;
    uint32_t n = local_args.N;
    uint32_t k = local_args.K;
    precision_t prec = (precision_t)local_args.prec;
    uint32_t setup_ssr = local_args.setup_ssr;
    uint32_t m_tiles = local_args.m_tiles;
    uint32_t n_tiles = local_args.n_tiles;
    uint32_t transa = local_args.transa;
    uint32_t transb = local_args.transb;
    double alpha = local_args.alpha;
    void* a = local_args.a;
    void* b = local_args.b;
    uint32_t beta = local_args.beta;
    void* c = local_args.c;

    // Calculate tile sizes
    uint32_t frac_m = m / m_tiles;
    uint32_t frac_n = n / n_tiles;
    uint32_t frac_k = k;
    uint32_t frac_a = frac_m * frac_k;
    uint32_t frac_c = frac_m * frac_n;
    uint32_t size_frac_a = frac_a * prec;
    uint32_t size_frac_b = frac_k * frac_n * prec;
    uint32_t size_frac_c = frac_c * prec;

    // Allocate A, B and C (double) buffers in TCDM
    void *local_a[2];
    void *local_b[2];
    void *local_c[2];

    // void* heap_ptr = (void*)(ALIGN_UP_TCDM((int)&local_args + sizeof(gemm_args_t)));
    void* heap_ptr = snrt_l1_next();
    int banks_per_buffer = snrt_cluster_compute_core_num();

    // The A, B and C buffers are stored in separate banks.
    // Particularly, every buffer spans a contiguous set of banks,
    // as many as the number of compute cores, which would access
    // them in parallel.
    local_a[0] = heap_ptr;
    local_b[0] = local_a[0] + BANK_WIDTH * banks_per_buffer;
    local_c[0] = local_b[0] + BANK_WIDTH * banks_per_buffer;
    if (SNRT_TCDM_HYPERBANK_NUM == 2) {
        local_a[1] = local_a[0] + SNRT_TCDM_HYPERBANK_SIZE;
        local_b[1] = local_b[0] + SNRT_TCDM_HYPERBANK_SIZE;
        local_c[1] = local_c[0] + SNRT_TCDM_HYPERBANK_SIZE;
    } else if (SNRT_TCDM_BANK_NUM == 64) {
        local_a[1] = local_c[0] + BANK_WIDTH * banks_per_buffer;
        local_b[1] = local_a[1] + BANK_WIDTH * banks_per_buffer;
        local_c[1] = local_b[1] + BANK_WIDTH * banks_per_buffer;
    } else {
        local_a[1] = local_a[0] + SNRT_TCDM_SIZE / 2;
        local_b[1] = local_b[0] + SNRT_TCDM_SIZE / 2;
        local_c[1] = local_c[0] + SNRT_TCDM_SIZE / 2;
    }

    // Calculate number of iterations
    int iterations = m_tiles * n_tiles + 2;
    int buff_idx;
    int i, i_dma_out, i_dma_in, i_compute;

    // Iterate over all tiles
    for (i = 0; i < iterations; i++) {
        if (snrt_is_dm_core()) {
            // DMA out
            // (out before in to avoid overwriting data)
            if (i > 1) {

                // Compute tile and buffer indices
                i_dma_out = i - 2;
                buff_idx = i_dma_out % 2;

                // Copy job outputs from TCDM
                snrt_dma_start_2d_wideptr(
                    c + i_dma_out * size_frac_c,
                    local_c[buff_idx],
                    BANK_WIDTH * banks_per_buffer,
                    BANK_WIDTH * banks_per_buffer,
                    HYPERBANK_WIDTH,
                    size_frac_c / (BANK_WIDTH * banks_per_buffer)
                );
                snrt_dma_wait_all();
            }

            // DMA in
            if (i < m_tiles * n_tiles) {

                // Compute tile and buffer indices
                i_dma_in = i;
                buff_idx = i_dma_in % 2;

                // Copy job operands in TCDM
                snrt_dma_start_2d_wideptr(
                    local_a[buff_idx],
                    a + i_dma_in * size_frac_a,
                    BANK_WIDTH * banks_per_buffer,
                    HYPERBANK_WIDTH,
                    BANK_WIDTH * banks_per_buffer,
                    size_frac_a / (BANK_WIDTH * banks_per_buffer)
                );
                snrt_dma_start_2d_wideptr(
                    local_b[buff_idx],
                    b,
                    BANK_WIDTH * banks_per_buffer,
                    HYPERBANK_WIDTH,
                    BANK_WIDTH * banks_per_buffer,
                    size_frac_b / (BANK_WIDTH * banks_per_buffer)
                );
                snrt_dma_wait_all();
            }
        }

        // Compute
        if (snrt_is_compute_core()) {
            if (i > 0 && i < (m_tiles * n_tiles + 1)) {

                // Compute tile and buffer indices
                i_compute = i - 1;
                buff_idx = i_compute % 2;

                // Perform tile computation
                volatile uint32_t ldb = frac_n;
                volatile uint32_t ldc = frac_n;
                if (transb) {
                    ldb = frac_k;
                }
                if (i_compute != 0) local_args.setup_ssr = 0;
                snrt_mcycle();
                sc_st_gemm(&local_args, local_a[buff_idx], local_b[buff_idx], beta,
                           local_c[buff_idx]);
                snrt_mcycle();
            }
        }

        // Synchronize cores after first iteration, exclusively for benchmarking
        snrt_cluster_hw_barrier();
    }

    return 0;
}