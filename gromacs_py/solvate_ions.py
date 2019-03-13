#!/usr/bin/env python3

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Solvate a gromacs system with water and add ions to neutralize the system charge and to reach an ionic concentration")
    parser.add_argument('-f', action="store", dest="f", help='Input PDB file', type=str, required=True)
    parser.add_argument('-p', action="store", dest="p", help='Topologie in gromacs format .top', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o", help='Output Directory', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name", help='Output file name', type=str, required=True)
    parser.add_argument('-d', action="store", dest="dist", help='Distance between the solute and the box', type=float, default=1.1)
    parser.add_argument('-C', action="store", dest="Conc", help='Ion concentration (mM), default = 0.15 (150mM)', type=float, default=0.15)
    return(parser)


if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()

    sys_top = gmx.GmxSys(name=args.name, coor_file=args.f, top_file=args.p)

    sys_top.solvate_add_ions(out_folder=args.o, name=args.name, ion_C=args.Conc, box_dist=args.dist)

    print("\n\nTopologie creation was sucessfull \n\tTopologie directorie :\t" + args.o)
    sys_top.display()
