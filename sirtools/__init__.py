"""
sirtools package.

This package provides utilities for fitting selective inversion relaxation
experiment data using CIFIT-style mechanism and data files.
"""

from .sirfit import do_fit, main, read_dat_file, read_mch_file

__all__ = ["main", "read_mch_file", "read_dat_file", "do_fit"]
