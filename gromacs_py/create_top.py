#!/usr/bin/env python3

""" Create topologie from a pdb file
"""

import argparse
import gromacs.gmx5 as gmx

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Create the topologie file\
     from a structure pdb file")
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o",
                        help='Output directory', type=str, required=True)
    parser.add_argument('-vsite', action="store_true", dest="vsite_flag",
                        help='Use virtual site for hydrogens')

    return parser


if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()
    if args.vsite_flag:
        vsite = "hydrogens"
    else:
        vsite = "none"

    sys_name = args.f.split("/")[-1][:-4]

    md_sys = gmx.GmxSys(name=sys_name, coor_file=args.f)
    md_sys.prepare_top(out_folder=args.o, vsite=vsite)
    md_sys.create_box(dist=1.0, box_type="dodecahedron", check_file_out=True)

    print("\n\nTopologie creation was sucessfull \n\tTopologie directorie :\t" + args.o)

    md_sys.display()
