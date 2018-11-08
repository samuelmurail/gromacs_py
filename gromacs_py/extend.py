#!/usr/bin/env python3

# Minimize a pdb file, return a pdb and a topologie

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse
from glob import glob
import os

def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Extend Simulation production/equilibration")
    parser.add_argument('-s', action="store", dest="tpr", help='Input tpr', type=str, required=True)
    parser.add_argument('-time', action="store", dest="time", help='Extend simulation time, default=10', type=float, default=10)
    parser.add_argument('-dt', action="store", dest="dt", help='integration time step, default=0.005', type=float, default=0.005)

    parser.add_argument('-nt', action="store", dest="nt", help='Total number of threads to start, default=0', type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi", help='Number of thread-MPI threads to start, default=0', type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid", help='List of GPU device id-s to use, default=\"\" ', default="None")
    return(parser)

if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()
    
    dt = args.dt
    nsteps = 1000*args.time/dt
    
    sys_prod = gmx.Gmx_sys()
    sys_prod.nt = args.nt
    sys_prod.ntmpi = args.ntmpi
    sys_prod.gpuid = args.gpuid
    
    sys_prod.extend_equi_prod(args.tpr, nsteps = nsteps)
    

    print("\n\nProduction extension was sucessfull ")
    
    sys_prod.display()
