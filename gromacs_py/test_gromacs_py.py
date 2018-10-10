#!/usr/bin/env python3

import doctest
import shutil
import os 

# Module to test:
import tools.pdb_manip
import tools.pdb2pqr
import gromacs.gmx5


test_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(test_path+"/tools")
shutil.rmtree('../test/output/pdb_manip_test', ignore_errors=True)
shutil.rmtree('../test/output/pdb2pqr_test', ignore_errors=True)

print("tools.pdb_manip:\t",doctest.testmod(tools.pdb_manip))
print("tools.pdb2pqr:  \t",doctest.testmod(tools.pdb2pqr))

os.chdir("../gromacs")

shutil.rmtree('../test/output/gmx5', ignore_errors=True)

print("gromacs.gmx5:    \t",doctest.testmod(gromacs.gmx5))



