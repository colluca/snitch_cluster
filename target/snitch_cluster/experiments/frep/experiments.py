#!/usr/bin/env python3
# Copyright 2025 ETH Zurich and University of Bologna.
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0
#
# Luca Colagrande <colluca@iis.ee.ethz.ch>

from copy import deepcopy
from snitch.target.SimResults import SimRegion
from snitch.target.experiment_utils import ExperimentManager
import random
import pandas as pd

from mako.template import Template
from pathlib import Path

NUM_GEMM_SIZES = 50
HW_CFGS = [
    'base32fc',
    'zonl32fc',
    'zonl64fc',
    'zonl64dobu',
    'zonl48dobu',
]
NUM_TILES = 3
ROI = 'tile_1'

POWER_GROUPS = [
    '*i_snitch_shared_muldiv',
    # '*i_snitch_cc',
    '*i_snitch_fp_ss*',
    '*i_idma_inst64*',
    '*i_snitch_icache*',
    '*i_data_mem*',
    '*i_axi_dma_xbar',
    '*i_axi_zeromem',
]

AREA_GROUPS = {
    'muldiv': '*i_snitch_shared_muldiv',
    # '*i_snitch_cc',
    'fpu': '*i_snitch_fp_ss*',
    'dma': '*i_idma_inst64*',
    'icache': '*i_snitch_icache*',
    'tcdm': '*i_data_mem*',
}


VSIM_BINS = {
    'base32fc':   str(Path.cwd() / 'hw/base32fc/bin/snitch_cluster.vsim'),
    'zonl32fc':   str(Path.cwd() / 'hw/zonl32fc/bin/snitch_cluster.vsim'),
    'zonl64fc':   str(Path.cwd() / 'hw/zonl64fc/bin/snitch_cluster.vsim'),
    'zonl64dobu': str(Path.cwd() / 'hw/zonl64dobu/bin/snitch_cluster.vsim'),
    'zonl48dobu': str(Path.cwd() / 'hw/zonl48dobu/bin/snitch_cluster.vsim'),
}
DATA_DIR = Path('data').absolute()
VERIFY_PY = Path('../../../../sw/blas/gemm/scripts/verify.py').absolute()


class FrepExperimentManager(ExperimentManager):

    def derive_axes(self, experiment):
        return {
            'hw': experiment['hw'],
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

    def derive_hw_cfg(self, experiment):
        return Path.cwd() / 'cfg' / f'{experiment["hw"]}.hjson'

    def derive_cdefines(self, experiment):
        if experiment['hw'] == 'base32fc':
            return {'USE_NESTED_FREP': 0}
        else:
            return {}


def generate_mat_size():
    sizes = list(range(8, 257, 8))
    return random.choice(sizes)


def calculate_total_size(m, n, k):
    prec = 8
    a_size = m * k * prec
    b_size = k * n * prec
    c_size = m * n * prec
    return 2 * (a_size + b_size + c_size)


def gen_experiments():
    # new layout requires checking that individual matrices fit in 8banks
    BANK_SIZE = 1 * 1024  # 2KiB
    MAX_ALLOWED_SIZE = 8 * BANK_SIZE  # Every matrix can take up maximum 8 banks

    experiments = []
    # TODO: filter repeated experiments
    while len(experiments) < NUM_GEMM_SIZES:
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

    return experiments


def get_average_fpu_util(row):
    util = []
    for i in range(0, 8):
        util.append(row['results'].get_metric(SimRegion(f'hart_{i}', ROI), 'fpss_fpu_occupancy'))
    return sum(util) / len(util)


def get_total_power(row):
    return row['power_results'].total_power


def get_area(row, key):
    return row['area_results'].cell_area[key]

def get_hier_area(row, val):
    df = row['area_groups']
    print(df[df['Name'] == val]['Area'].item())
    return df[df['Name'] == val]['Area'].item()

def main():

    seed = 31

    random.seed(seed)

    # m, n and k are the dimensions of the tile
    experiments = gen_experiments()
    experiments = [deepcopy(exp) for _ in range(len(HW_CFGS)) for exp in experiments]

    for i in range(len(HW_CFGS)):
        for j in range(NUM_GEMM_SIZES):
            k = i * NUM_GEMM_SIZES + j
            experiments[k]['hw'] = HW_CFGS[i]
            experiments[k]['app'] = 'gemm'
            experiments[k]['m_tiles'] = NUM_TILES
            experiments[k]['n_tiles'] = 1
            experiments[k]['cmd'] = [str(VERIFY_PY), VSIM_BINS[HW_CFGS[i]], "${elf}"]

    manager = FrepExperimentManager(experiments)
    manager.run()
    # TODO: revisit the start and end times, probably makes sense to use the cluster barrier
    # between tiles as a delimiter
    manager.export_power_experiments(SimRegion('hart_0', ROI))

    df = manager.get_results()
    df['size (KiB)'] = df.apply(lambda row: calculate_total_size(row['m'], row['n'], row['k']) / 1024, axis=1)

    if manager.perf_results_available:
        df['fpu_util'] = df.apply(get_average_fpu_util, axis=1)
    # print(df)

    if manager.power_results_available:
        df['total_power'] = df.apply(get_total_power, axis=1)
        # print(df['total_power'])
        df.to_csv('power.csv', index=False)
        # breakdown = df.iloc[0]['power_results'].group_power_breakdown(POWER_GROUPS)
        # print(breakdown)
        # print(df.iloc[0]['power_results'].total_power)
        # sum = breakdown['total_power'].sum() + df.iloc[0]['power_results'].clock_power
        # print(
        #     breakdown['total_power'].sum(), '+',
        #     df.iloc[0]['power_results'].clock_power, '=',
        #     sum
        # )

    df_area = manager.get_area_results()
    if manager.area_results_available:
        for key in ['TotArea', 'CombArea', 'SeqArea', 'MacroArea', 'BufInvArea']:
            df_area[key] = df_area.apply(lambda row: get_area(row, key), axis=1)

        df_area['area_groups'] = df_area.apply(lambda row: row['area_results'].group_area_breakdown(AREA_GROUPS.values()), axis=1)
        for key, val in AREA_GROUPS.items():
            print(key)
            df_area[key] = df_area.apply(lambda row: get_hier_area(row, val), axis=1)

        df_area['area_groups'] = None
        print(df_area)
        df_area.to_csv('area.csv', index=False)

    # Export results to file
    df.to_csv('results.csv', index=False)


if __name__ == '__main__':
    main()
