#!/usr/bin/env python3

# Minimize a pdb file, return a pdb and a topologie

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Minimize a pdb structure")
    parser.add_argument('-fsys', action="store", dest="f_sys", help='Input PDB file of the system')
    parser.add_argument('-psys', action="store", dest="p_sys", help='Topologie in gromacs format .top of the system')
    parser.add_argument('-fmol', action="store", dest="f_mol", help='Input PDB file of the molecule to insert')
    parser.add_argument('-pmol', action="store", dest="p_mol", help='Topologie in gromacs format .top of the molecule to insert')
    parser.add_argument('-nmol', action="store", dest="num_mol", help='Number of molecule to insert', type=int, default=20)
    parser.add_argument('-o', action="store", dest="o", help='Output Directory')
    parser.add_argument('-n', action="store", dest="name", help='Output file name')

    return parser


if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()

    sys_raw = gmx.GmxSys(name=args.name, coor_file=args.f_sys, top_file=args.p_sys)
    mol_gmx = gmx.GmxSys(name="mol", coor_file=args.f_mol, top_file=args.p_mol)

    sys_raw.display()
    mol_gmx.display()

    sys_raw.insert_mol_sys(mol_gromacs=mol_gmx, mol_num=args.num_mol, new_name=args.name, out_folder=args.o, check_file_out=True)

    print("\n\nInsertion was sucessfull \n\tSystem directory :\t" + args.o + "\n\tsystem coor file:\t" + sys_raw.coor_file + "\n\tsystem top file:\t" + sys_raw.top_file)
