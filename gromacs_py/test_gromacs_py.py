#!/usr/bin/env python3

""" Script to test core functions and classes of gromacs_py
"""

import doctest
import shutil
import sys
import os

# Module to test:
import gromacs.tools.pdb_manip
import gromacs.tools.os_command
import gromacs.tools.pdb2pqr
import gromacs.gmx5

__author__ = "Samuel Murail"


print("gromacs location package :",gromacs.gmx5.__file__)
print("tools.os_command:  \t", doctest.testmod(gromacs.tools.os_command))
print("tools.pdb_manip:\t", doctest.testmod(gromacs.tools.pdb_manip))
print("tools.pdb2pqr:  \t", doctest.testmod(gromacs.tools.pdb2pqr))
print("gromacs.gmx5:    \t", doctest.testmod(gromacs.gmx5))

shutil.rmtree('gromacs_py_test_out', ignore_errors=True)
