#!/usr/bin/env python3
# Copyright 2025 ETH Zurich and University of Bologna.
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0
#
# Luca Colagrande <colluca@iis.ee.ethz.ch>

from snitch.target.SimResults import SimRegion
from snitch.target.experiment_utils import ExperimentManager
import random

from mako.template import Template
from pathlib import Path

DATA_DIR = Path('data').absolute()
VERIFY_PY = Path('../../../../sw/blas/gemm/scripts/verify.py').absolute()

class FrepExperimentManager(ExperimentManager):

    def derive_axes(self, experiment):
        return {
            'm': experiment['m'],
            'n': experiment['n'],
            'k': experiment['k'],
        }

    def derive_data_cfg(self, experiment):
        # Create parent directory for configuration file
        cfg_path = DATA_DIR / experiment['name'] / 'cfg.json'
        cfg_path.parent.mkdir(parents=True, exist_ok=True)

        # Fill in configuration template and write configuration file
        with open('cfg.json.tpl') as f:
            cfg = Template(f.read()).render(experiment=experiment)
        with open(cfg_path, 'w') as f:
            f.write(cfg)
        return cfg_path


def generate_mat_size():
    sizes = list(range(8, 257, 8))
    return random.choice(sizes)


def calculate_total_size(m, n, k):
    prec = 8
    a_size = m * k * prec
    b_size = k * n * prec
    c_size = m * n * prec
    return 2 * (a_size + b_size + c_size)


def gen_experiment(experiments):
    # new layout requires checking that individual matrices fit in 8banks
    BANK_SIZE = 1 * 1024  # 2KiB
    MAX_ALLOWED_SIZE = 8 * BANK_SIZE  # Every matrix can take up maximum 8 banks

    i = 0
    while i < 10:
        # generate random values in [8, 16, 24, ..., 256] for m, n, k
        m = generate_mat_size()
        n = generate_mat_size()
        k = generate_mat_size()

        # check if the matrices fit in TCDM
        # new layout requires checking that individual matrices fit in 8banks
        prec = 8
        a_size = m * k * prec
        b_size = k * n * prec
        c_size = m * n * prec
        max_size = max(a_size, b_size, c_size)

        if max_size < MAX_ALLOWED_SIZE:
            experiments.append({'m': m, 'n': n, 'k': k})
            i += 1

    return experiments


def main():

    seed = 31

    random.seed(seed)

    # m, n and k are the dimensions of the tile
    experiments = []

    experiments = gen_experiment(experiments)

    for experiment in experiments:
        experiment['app'] = 'gemm'
        experiment['m_tiles'] = 2
        experiment['n_tiles'] = 1
        experiment['cmd'] = [str(VERIFY_PY), "${sim_bin}", "${elf}"]

    manager = FrepExperimentManager(experiments)
    manager.run()
    manager.export_power_experiments(SimRegion('hart_0', 'tile_1'))

    df = manager.get_results()
    df['fpu_util'] = df.apply(lambda row: row['results'].get_metric(SimRegion('hart_0', 'tile_1'), 'fpss_fpu_occupancy'), axis=1)
    df['size (KiB)'] = df.apply(lambda row: calculate_total_size(row['m'], row['n'], row['k']) / 1024, axis=1)
    print(df)


if __name__ == '__main__':
    main()
