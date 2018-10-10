#!/usr/bin/env python3

##################################
#########   TEST PART   ##########
##################################

__author__ = "Samuel Murail"

import sys
import os

# necessary to load modules in root dir 
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import gromacs.gmx5 as gmx
import tools.pdb_manip as pdb
import tools.osCommand as osCommand



if __name__ == '__main__':


    dir_out = "./output/"
    
    

    CYPA = pdb.coor()
    CYPA.read_pdb("input/1AWR-055_bestene1-mc.pdb")
    CYPA.write_pdb(dir_out+"tmp.pdb")