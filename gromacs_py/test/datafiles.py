#!/usr/bin/env python3
# coding: utf-8

"""
Location of datafiles for gromacs_py tests

Use it as:

```
from gromacs_py.test.datafiles import *
```
"""

import os

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

PYTEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_FILE_PATH = os.path.join(PYTEST_DIR, "../test_files/")

PDB_1D30 = os.path.join(TEST_FILE_PATH, '1D30.pdb')
PDB_1Y0M = os.path.join(TEST_FILE_PATH, '1y0m.pdb')
PDB_1RXZ = os.path.join(TEST_FILE_PATH, '1rxz.pdb')
PDB_1JD4 = os.path.join(TEST_FILE_PATH, '1jd4.pdb')
