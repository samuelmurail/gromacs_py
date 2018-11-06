#!/usr/bin/env python3

import doctest
import shutil
import os 

# Module to test:
import tools.pdb_manip
import tools.os_command
import tools.pdb2pqr
import gromacs.gmx5



print("tools.os_command:  \t",doctest.testmod(tools.os_command))
print("tools.pdb_manip:\t",doctest.testmod(tools.pdb_manip))
print("tools.pdb2pqr:  \t",doctest.testmod(tools.pdb2pqr))
print("gromacs.gmx5:    \t",doctest.testmod(gromacs.gmx5))

shutil.rmtree('gromacs_py_test_out', ignore_errors=True)




