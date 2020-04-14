#!/usr/bin/env python3

""" Script to test core functions and classes of gromacs_py
"""

import doctest
import shutil

# Module to test:
import gromacs.gmx5

__author__ = "Samuel Murail"


print("gromacs location package :",gromacs.gmx5.__file__)
print("gromacs.gmx5:    \t", doctest.testmod(gromacs.gmx5))

shutil.rmtree('gromacs_py_test_out', ignore_errors=True)
