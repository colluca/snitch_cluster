# Copyright 2023 ETH Zurich and University of Bologna.
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0
#
# Author: Luca Colagrande <colluca@iis.ee.ethz.ch>
"""Convenience functions and classes for data generation scripts."""


import argparse
import json5
import pathlib
import struct
from datetime import datetime
import torch
import numpy as np
import pyflexfloat as ff
import humanize
import math

# Maximum available size in TCDM (in bytes)
TCDM_HEAP_SIZE = 112 * 1024


def emit_license():
    """Emit license header.

    Returns:
        A header string.
    """

    s = (f"// Copyright {datetime.now().year} ETH Zurich and University of Bologna.\n"
         f"// Licensed under the Apache License, Version 2.0, see LICENSE for details.\n"
         f"// SPDX-License-Identifier: Apache-2.0\n")
    return s


# Enum value can be a string or an integer, this function uniformizes the result to integers only
def _integer_precision_t(prec):
    if isinstance(prec, str):
        return {'FP64': 8, 'FP32': 4, 'FP16': 2, 'FP8': 1}[prec]
    else:
        return prec


def size_from_precision_t(prec):
    """Return the size in bytes of a `precision_t` type.

    Args:
        prec: A value of type `precision_t`. Accepts both enum strings
            (e.g. "FP64") and integer enumeration values (e.g. 8).
    """
    return _integer_precision_t(prec)


def ff_desc_from_precision_t(prec):
    """Convert `precision_t` type to a FlexFloat descriptor.

    Args:
        prec: A value of type `precision_t`. Accepts both enum strings
            (e.g. "FP64") and integer enumeration values (e.g. 8).
    """
    precision_t_to_ff_desc_map = {
        8: 'fp64',
        4: 'fp32',
        2: 'fp16',
        1: 'e5m2'
    }
    return precision_t_to_ff_desc_map[_integer_precision_t(prec)]


def torch_type_from_precision_t(prec):
    """Convert `precision_t` type to PyTorch type.

    Args:
        prec: A value of type `precision_t`. Accepts both enum strings
            (e.g. "FP64") and integer enumeration values (e.g. 8).
    """
    precision_t_to_torch_type_map = {
        8: torch.float64,
        4: torch.float32,
        2: torch.float16,
        1: torch.float8_e4m3fn
    }
    return precision_t_to_torch_type_map[_integer_precision_t(prec)]


# Returns the C type representing a floating-point value of the specified precision
def ctype_from_precision_t(prec):
    """Convert `precision_t` type to a C type string.

    Args:
        prec: A value of type `precision_t`. Accepts both enum strings
            (e.g. "FP64") and integer enumeration values (e.g. 8).
    """
    precision_t_to_ctype_map = {
        8: 'double',
        4: 'float',
        2: '__fp16',
        1: '__fp8'
    }
    return precision_t_to_ctype_map[_integer_precision_t(prec)]


def flatten(array):
    """Flatten various array types with a homogeneous API.

    Args:
        array: Can be a Numpy array, a PyTorch tensor or a nested list.
    """
    if isinstance(array, np.ndarray):
        return array.flatten()
    elif isinstance(array, torch.Tensor):
        return array.numpy().flatten()
    elif isinstance(array, list):
        return np.array(array).flatten()


def _variable_attributes(alignment=None, section=None):
    attributes = ''
    if alignment:
        attributes = f'__attribute__ ((aligned ({alignment})))'
    if section:
        attributes += f' __attribute__ ((section ("{section}")))'
    return attributes


def _alias_dtype(dtype):
    if dtype == '__fp8':
        return 'char'
    else:
        return dtype


def format_array_declaration(dtype, uid, shape, alignment=None, section=None):
    attributes = _variable_attributes(alignment, section)
    s = f'{_alias_dtype(dtype)} {uid}'
    for dim in shape:
        s += f'[{dim}]'
    if attributes:
        s += f' {attributes};'
    else:
        s += ';'
    return s


# In the case of dtype __fp8, array field expects a dictionary of
# sign, exponent and mantissa arrays
def format_array_definition(dtype, uid, array, alignment=None, section=None):
    # Definition starts with the declaration stripped off of the terminating semicolon
    s = format_array_declaration(dtype, uid, array.shape, alignment, section)[:-1]
    s += ' = '
    s += format_array_initializer(dtype, array)
    s += ';'
    return s


def format_scalar_definition(dtype, uid, scalar):
    s = f'{_alias_dtype(dtype)} {uid} = {scalar};'
    return s


def format_array_initializer(dtype, array):
    s = '{\n'
    array = flatten(array)
    for el in array:
        if dtype == '__fp8':
            el_str = f'{hex(el.bits())}'
        else:
            el_str = f'{el}'
        s += f'\t{el_str},\n'
    s += '}'
    return s


def format_struct_definition(dtype, uid, map):
    def format_value(value):
        if isinstance(value, list):
            return format_array_initializer('str', value)
        elif isinstance(value, bool):
            return str(int(value))
        else:
            return str(value)

    filtered_map = {key: value for key, value in map.items() if value is not None and value != ''}

    formatted_items = [
        f'\t.{key} = {format_value(value)}'
        for key, value in filtered_map.items()
    ]
    s = f'{_alias_dtype(dtype)} {uid} = {{\n'
    s += ',\n'.join(formatted_items)
    s += '\n};'
    return s


def format_ifdef_wrapper(macro, body):
    s = f'#ifdef {macro}\n'
    s += f'{body}\n'
    s += f'#endif // {macro}\n'
    return s


def to_bits(chunks, ctype='uint32_t'):
    """Convert structured data to a bit array.

    Returns a list of bits corresponding to the structured data.
    """
    if ctype == 'uint32_t':
        bitstring = 0

        # Combine chunks into a single bitstring
        for i, chunk in enumerate(chunks):
            bitstring |= int(chunk) << (32 * i)

        # Extract individual bits from the bitstring
        bit_array = [(bitstring >> i) & 1 for i in range(bitstring.bit_length())]
        return bit_array
    else:
        raise ValueError(f'Unsupported ctype: {ctype}')


def from_bits(bit_array, ctype='uint32_t'):
    """Get structured data from a list of bits.

    Returns a homogeneous list of the specified type from the list of
    bits.
    """
    if ctype == 'uint32_t':
        # Build bitstring
        bitstring = 0
        for i, bit in enumerate(bit_array):
            bitstring |= bit << i

        # Split bitstring into chunks and convert to specified type
        num_chunks = math.ceil(math.log2(bitstring) / 32)
        chunks = [(bitstring >> (32 * i)) & 0xFFFFFFFF for i in range(num_chunks)]
        return chunks
    else:
        raise ValueError(f'Unsupported ctype: {ctype}')


def from_buffer(byte_array, ctype='uint32_t'):
    """Get structured data from raw bytes.

    If `ctype` is a C type string, it returns a homogeneous list of the
    specified type from the raw data, using numpy's `from_buffer` method.
    Note that if `ctype` is equal to `__fp8`, given that there is no
    native fp8 format in Numpy, an array of FlexFloat objects is returned.

    Alternatively, a dictionary can be passed to `ctype` to extract a
    struct from the raw data. In this case, it returns a dictionary with
    the same keys as in `ctype`. The values in the `ctype` dictionary
    should be format strings compatible with Python's `struct` library.
    The order of the keys in the `ctype` dictionary should reflect the
    order in which the variables appear in the raw data.
    """
    # Types which have a direct correspondence in Numpy
    NP_DTYPE_FROM_CTYPE = {
        'uint32_t': np.uint32,
        'int32_t': np.int32,
        'int': np.int32,
        'double': np.float64,
        'float': np.float32,
        '__fp16': np.float16
    }

    if isinstance(ctype, dict):
        # byte_array assumed little-endian
        struct_fields = ctype.keys()
        fmt_specifiers = ctype.values()
        fmt_string = ''.join(fmt_specifiers)
        field_values = struct.unpack(f'<{fmt_string}', byte_array)
        return dict(zip(struct_fields, field_values))
    elif ctype in NP_DTYPE_FROM_CTYPE.keys():
        dtype = NP_DTYPE_FROM_CTYPE[ctype]
        return np.frombuffer(byte_array, dtype=dtype)
    elif ctype == '__fp8':
        return ff.frombuffer(byte_array, 'e5m2')


class DataGen:
    """Base data generator class.

    Base data generator class which can be inherited to easily develop a
    custom data generator script for any kernel.
    """

    def parser(self):
        """Default argument parser for data generation scripts.

        It is an instance of the `ArgumentParser` class from the `argparse`
        module. Subclasses can extend this and add custom arguments via
        the parser's `add_argument()` method.
        """
        parser = argparse.ArgumentParser(description='Generate data for kernels')
        parser.add_argument(
            "-c", "--cfg",
            type=pathlib.Path,
            required=True,
            help='Select param config file kernel'
        )
        parser.add_argument(
            '--section',
            type=str,
            help='Section to store matrices in')
        parser.add_argument(
            'output',
            type=pathlib.Path,
            help='Path of the output header file')
        return parser

    def parse_args(self):
        """Parse default data generation script arguments.

        Returns the arguments passed to the data generation script, parsed
        using the `parser()` method.
        """
        return self.parser().parse_args()

    def emit_header(self, **kwargs):
        """Emits a C header containing generated data.

        The base implementation emits a string which only contains a
        license header. Subclasses should extend this method and append
        the generated data to the license header.

        Returns:
            A string with the generated C header contents.
        """
        return emit_license()

    def main(self):
        """Default main function for data generation scripts."""
        self.args = self.parse_args()

        # Load param config file
        with self.args.cfg.open() as f:
            param = json5.loads(f.read())
        param['section'] = self.args.section

        # Emit header file
        with open(self.args.output, 'w') as f:
            f.write(self.emit_header(**param))


def validate_tcdm_footprint(size, silent=False):
    """Check whether data of specified size fits in TCDM.

    Throws an assertion error if the specified size exceeds the space
    available for the heap in TCDM.

    Args:
        size: The size of the data in bytes.
        silent: If True, will not print the size to stdout.
    """
    assert size < TCDM_HEAP_SIZE, \
        f'Total heap space required {humanize.naturalsize(size, binary=True)} exceeds ' \
        f'limit of {humanize.naturalsize(TCDM_HEAP_SIZE, binary=True)}'
    if not silent:
        print(f'Total heap space required {humanize.naturalsize(size, binary=True)}')
