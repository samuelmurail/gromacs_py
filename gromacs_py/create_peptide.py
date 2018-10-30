#!/usr/bin/env python3

# Create a Peptide

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse
from glob import glob
import os


def parser_input():

    parser = argparse.ArgumentParser(description="Create a linear peptide structure, do a minimisation and a vacuum equilibration")
    parser.add_argument('-seq', action="store", dest="seq", help='Peptide sequence', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o", help='Output Directory', type=str, required=True)
    parser.add_argument('-m_steps', action="store", dest="min_steps", help='Minimisation nsteps, default=1000', type=int, default=1000)
    parser.add_argument('-time', action="store", dest="time", help='Vacuum equilibration time(ns), default = 1ns', type=float, default=1)

    return(parser)




if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()

    dt = 0.001
    step = 1000*args.time/dt
    sequence = args.seq
    out_folder = args.o
    em_nsteps = args.min_steps
    
    peptide =  gmx.gmx_sys(name='pep_'+sequence)
    peptide.create_peptide(sequence = sequence, out_folder = out_folder, em_nsteps = em_nsteps, equi_nsteps = step)
    
    print("\n\nPeptide Creation was sucessfull \n\tPeptide directorie :\t"+args.o)
    peptide.display()
