#!/usr/bin/env python3

""" Script to test core functions and classes of gromacs_py
"""

import doctest
import shutil
import sys
import os

# In order to fix issue with relative import like ..tools
# the test has to load gromacs_py.gromacs.gmx5 and not gromacs.gmx5
# (Can't load a module outside its package)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import gromacs_py.tools.pdb_manip
import gromacs_py.tools.os_command
import gromacs_py.tools.pdb2pqr
import gromacs_py.gromacs.gmx5


__author__ = "Samuel Murail"


#sys.path.append("..")
print("gromacs location package :",gromacs_py.__file__)
print("tools.os_command:  \t", doctest.testmod(gromacs_py.tools.os_command))
print("tools.pdb_manip:\t", doctest.testmod(gromacs_py.tools.pdb_manip))
print("tools.pdb2pqr:  \t", doctest.testmod(gromacs_py.tools.pdb2pqr))
print("gromacs.gmx5:    \t", doctest.testmod(gromacs_py.gromacs.gmx5))

shutil.rmtree('gromacs_py_test_out', ignore_errors=True)
