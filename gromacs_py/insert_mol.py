#!/usr/bin/env python3

# Minimize a pdb file, return a pdb and a topologie

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse
from glob import glob
import os

# Parse arguments :
parser = argparse.ArgumentParser(description="Minimize a pdb structure")
parser.add_argument('-fsys', action="store", dest="f_sys", help='Input PDB file of the system')
parser.add_argument('-psys', action="store", dest="p_sys", help='Topologie in gromacs format .top of the system')
parser.add_argument('-fmol', action="store", dest="f_mol", help='Input PDB file of the molecule to insert')
parser.add_argument('-pmol', action="store", dest="p_mol", help='Topologie in gromacs format .top of the molecule to insert')
parser.add_argument('-nmol', action="store", dest="num_mol", help='Number of molecule to insert', type=int, default=20)
parser.add_argument('-o', action="store", dest="o", help='Output Directory')
parser.add_argument('-n', action="store", dest="name", help='Output file name')
args = parser.parse_args()




#print("Min steps :\t",args.min_steps,"\nEqui HA time :",args.HA_time,"ns\nEqui CA time :",args.CA_time,"ns\nEqui CA_LOW time :",args.CA_LOW_time,"ns")

sys_min = gmx.insert_mol_sys_vmd(sys_pdb = args.f_sys, sys_top = args.p_sys, mol_pdb = args.f_mol, mol_top = args.p_mol,
	mol_num = args.num_mol, mol_length = 3, out_folder = args.o, sys_name = args.name)

### TO DO ####
print("\n\nInsertion was sucessfull \n\tSystem directory :\t"+args.o+"\n\tsystem coor file:\t"+args.o+"/"+args.name+"_neutral.pdb\n\tsystem top file:\t"+args.o+"/"+args.name+".top")

#~/Documents/repository/Gromacs_Setup/insert_mol.py 
#-fsys mini_equi_3_step_SOL_ION_2cyh/03_equi_CA_LOW/equi_CA_LOW_2cyh_compact.gro 
#-psys solv_ion_2cyh/solv_ion_2cyh.top 
#-fmol ARP/01_equi_imp/equi_implicit_compact.gro 
#-pmol ARP/ARP_pdb2gmx.top 
#-nmol 20 -o 2CYH_20_ARP -n 2cyh_20_ARP