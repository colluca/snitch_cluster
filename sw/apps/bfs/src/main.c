// Copyright 2024 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Author: Luca Colagrande <colluca@iis.ee.ethz.ch>

#include <stdint.h>

#include "data.h"

#include "bfs.h"

int main() {
    bfs_args_t args = {
        .num_vertices = num_vertices,
        .load_graph = 1,
        .graph_offsets_addr = (uint64_t)offsets,
        .graph_adjacencies_addr = (uint64_t)adjacencies,
        .frontier_addr = (uint64_t)frontier,
        .dist_addr = (uint64_t)dist,
        .out_frontier_addr = (uint64_t)out_frontier,
        .out_dist_addr = (uint64_t)out_dist
    };
    bfs_job(&args);
    return 0;
}
