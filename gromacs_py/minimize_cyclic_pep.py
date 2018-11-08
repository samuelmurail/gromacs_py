#!/usr/bin/env python3

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse
from glob import glob
import os
import shutil

def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description='Minimize a cyclic peptide structure in 2 steps, the first step without bonds constraints and the second step with bonds constraints')
    parser.add_argument('-f', action="store", dest="f", help='Input PDB file', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name", help='Output file name', type=str, required=True)
    parser.add_argument('-dir', action="store", dest="out_dir", help='Output directory for intermediate files', type=str, default="./tmp")
    parser.add_argument('-m_steps', action="store", dest="min_steps", help='Minimisation nsteps, default=1000', type=int, default=1000)
    parser.add_argument('-keep', action="store_true", dest="keep_flag", help='Flag to keep temporary files (without flag output directory will be delete')

    return(parser)


if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()


    vsite = "none"
    
    cyclic_pep = gmx.Gmx_sys(name = args.name, coor_file = args.f)
    cyclic_pep.cyclic_peptide_top(out_folder = args.out_dir+'/top')
    
    cyclic_pep.em(out_folder = args.out_dir+'/em',name = "min_"+args.name,
    	nsteps = args.min_steps, posres = "", nstxout = 1000,create_box_flag = True,
    	constraints = "none")
    
    cyclic_pep.convert_trj( traj=False )
    
    # Get the minimised structure:
    shutil.copyfile(cyclic_pep.coor_file, args.name+'.pdb')
    
    # Keep or not the intermediate files:
    if not args.keep_flag:
    	shutil.rmtree(args.out_dir, ignore_errors=True)
    else:
    	cyclic_pep.display()
    
