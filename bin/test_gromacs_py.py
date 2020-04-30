#!/usr/bin/env python3

""" Script to test core functions and classes of gromacs_py
"""

import doctest
import shutil

# Module to test:
from gromacs_py import gmx

__author__ = "Samuel Murail"


print("gromacs location package :", gmx.__file__)
print("gromacs.gmx5:    \t", doctest.testmod(gmx))

shutil.rmtree('gromacs_py_test_out', ignore_errors=True)
