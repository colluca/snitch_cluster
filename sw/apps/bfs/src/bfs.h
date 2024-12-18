// Copyright 2024 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Author: Luca Colagrande <colluca@iis.ee.ethz.ch>

#include <stdint.h>
#include <math.h>

#include "args.h"
#include "snrt.h"

__thread uint32_t size_of_int = sizeof(uint32_t);
__thread uint32_t *local_graph_offsets;
__thread uint32_t *local_graph_adjacencies;
__thread uint32_t *local_frontier;

void frontier_job(uint32_t num_vertices, uint32_t *graph_offsets, uint32_t *graph_adjacencies, uint32_t* frontier, int32_t* dist, uint32_t* new_frontier, int32_t* new_dist) {

    // Distribute workload
    uint32_t num_vertices_per_core = num_vertices / snrt_cluster_compute_core_num();

    // Iterate vertices
    for (uint32_t i = 0; i < num_vertices_per_core; i++) {

        // Get vertex index
        uint32_t v_cluster = snrt_cluster_core_idx() * num_vertices_per_core + i;
        uint32_t v = snrt_cluster_idx() * num_vertices + v_cluster;

        // Nothing to do if vertex has already been visited
        if (dist[v] == -1) {

            // Check if any of its neighbors are in the frontier
            for (uint32_t j = graph_offsets[v_cluster]; j < graph_offsets[v_cluster+1]; j++) {

                // Get vertex ID of neighbor
                uint32_t n = graph_adjacencies[j - graph_offsets[0]];

                // Index frontier vector using vertex ID of neighbor
                uint32_t neighbor_in_frontier = (frontier[n / 32] >> (n % 32)) & 1;

                // If neighbor is in the frontier, add current vertex to new frontier
                // and update its distance
                if (neighbor_in_frontier) {
                    new_dist[v_cluster] = dist[n] + 1;
                    __atomic_or_fetch(
                        &(new_frontier[v_cluster / 32]),
                        1 << (v_cluster % 32),
                        __ATOMIC_RELAXED
                    );
                    break;
                }
            }
        }
    }
}

void bfs_job(void *args) {
    uint32_t *local_new_frontier;
    int32_t *local_dist, *local_new_dist;
    bfs_args_t *local_args;

#ifndef JOB_ARGS_PRELOADED
    // Allocate space for job arguments in TCDM
    local_args = (bfs_args_t *)snrt_l1_alloc_cluster_local(
        sizeof(bfs_args_t), sizeof(double));

    // Copy job arguments to TCDM
    if (snrt_is_dm_core()) {
        snrt_dma_start_1d(local_args, args, sizeof(bfs_args_t));
        snrt_dma_wait_all();
    }
    snrt_cluster_hw_barrier();
#else
    local_args = (bfs_args_t *)args;
#endif

    // Aliases
    uint32_t num_vertices = local_args->num_vertices;
    uint32_t *frontier = (uint32_t *)(local_args->frontier_addr);
    int32_t *dist = (int32_t *)(local_args->dist_addr);
    uint32_t *out_frontier = (uint32_t *)(local_args->out_frontier_addr);
    int32_t *out_dist = (int32_t *)(local_args->out_dist_addr);

    // Allocate local variables
    uint32_t num_vertices_frac = num_vertices / snrt_cluster_num();
    // Ceil removed for performance reasons
    // size_t size_frontier = ceil(num_vertices / 32.0) * size_of_int;
    // size_t size_frontier_frac = ceil(num_vertices_frac / 32.0) * size_of_int;
    size_t size_frontier = num_vertices / 32 * sizeof(uint32_t);
    size_t size_frontier_frac = num_vertices_frac / 32 * sizeof(uint32_t);;
    size_t size_dist = num_vertices * sizeof(uint32_t);
    size_t size_dist_frac = num_vertices_frac * sizeof(uint32_t);
    if (local_args->load_graph) {
        uint32_t *graph_offsets = (uint32_t *)(local_args->graph_offsets_addr);
        uint32_t *graph_adjacencies = (uint32_t *)(local_args->graph_adjacencies_addr);

        size_t size_offsets = (num_vertices_frac + 1) * sizeof(uint32_t);
        local_graph_offsets = snrt_l1_alloc_cluster_local(size_offsets, sizeof(uint32_t));

        // Load offsets
        if (snrt_is_dm_core()) {
            snrt_dma_start_1d(
                local_graph_offsets,
                graph_offsets + snrt_cluster_idx() * num_vertices_frac,
                size_offsets
            );
            snrt_dma_wait_all();
        }
        snrt_cluster_hw_barrier();

        // Allocate memory needed for graph adjacencies
        uint32_t start_edge = local_graph_offsets[0];
        uint32_t num_edges = local_graph_offsets[num_vertices_frac] - start_edge;
        uint32_t size_adjacencies = num_edges * sizeof(uint32_t);
        local_graph_adjacencies = snrt_l1_alloc_cluster_local(size_adjacencies, sizeof(uint32_t));

        // Load adjacencies
        if (snrt_is_dm_core()) {
            snrt_dma_start_1d(
                local_graph_adjacencies,
                graph_adjacencies + start_edge,
                size_adjacencies);
        }
    }
    else {
        snrt_l1_update_next_v2(local_frontier);
    }
    local_frontier = snrt_l1_alloc_cluster_local(size_frontier, sizeof(uint32_t));
    local_new_frontier = snrt_l1_alloc_cluster_local(size_frontier_frac, sizeof(uint32_t));
    local_dist = snrt_l1_alloc_cluster_local(size_dist, sizeof(uint32_t));
    local_new_dist = snrt_l1_alloc_cluster_local(size_dist_frac, sizeof(uint32_t));

    // Load frontier and distance vectors
    if (snrt_is_dm_core()) {
        snrt_dma_start_1d(local_frontier, frontier, size_frontier);
        snrt_dma_start_1d(local_dist, dist, size_dist);
        snrt_dma_start_1d(local_new_frontier, (void *)snrt_zero_memory_ptr(), size_frontier_frac);
        snrt_dma_load_1d_tile(
            local_new_dist,
            local_dist,
            snrt_cluster_idx(),
            num_vertices_frac,
            sizeof(int32_t)
        );
        snrt_dma_wait_all();
    }

    // Synchronize with DM core to wait for job operands
    snrt_mcycle();
    snrt_cluster_hw_barrier();

    // Perform frontier calculation
    if (snrt_is_compute_core()) {
        snrt_mcycle();
        frontier_job(num_vertices_frac, local_graph_offsets, local_graph_adjacencies, local_frontier, local_dist, local_new_frontier, local_new_dist);
        snrt_mcycle();
    };

    // Synchronize with DM core to wait for job results
    snrt_cluster_hw_barrier();
    snrt_mcycle();

    // Copy data out of TCDM
    if (snrt_is_dm_core()) {
        snrt_dma_store_1d_tile(
            out_frontier,
            local_new_frontier,
            snrt_cluster_idx(),
            // Ceil removed for performance reasons
            // ceil(num_vertices_frac / 32.0),
            num_vertices_frac / 32,
            sizeof(uint32_t)
        );
        snrt_dma_store_1d_tile(
            out_dist,
            local_new_dist,
            snrt_cluster_idx(),
            num_vertices_frac,
            sizeof(int32_t)
        );
        snrt_dma_wait_all();
        snrt_mcycle();
    }

    snrt_cluster_hw_barrier();

    // Free memory
#ifndef JOB_ARGS_PRELOADED
    snrt_l1_update_next_v2(local_args);
#else
    snrt_l1_update_next_v2(local_graph_offsets);
#endif
}
