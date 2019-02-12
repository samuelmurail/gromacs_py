#!/usr/bin/env python3

""" Minimize a pdb file, return a pdb
"""

import argparse
import gromacs.gmx5 as gmx

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Minimize a pdb structure in 2 steps,\
                                     the first step without bonds constraints and the \
                                     second step with")
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', required=True)
    parser.add_argument('-p', action="store", dest="p",
                        help='Topologie in gromacs format .top', required=True)
    parser.add_argument('-o', action="store", dest="o",
                        help='Output Directory', required=True)
    parser.add_argument('-n', action="store", dest="name",
                        help='Output file name', required=True)
    parser.add_argument('-m_steps', action="store", dest="min_steps",
                        help='Minimisation nsteps, default=1000', type=int, default=1000)
    parser.add_argument('-box', action="store", dest="box",
                        help='Create a box, default=False', type=bool, default=False)
    parser.add_argument('-nt', action="store", dest="nt",
                        help='Total number of threads to start, default=0', type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi",
                        help='Number of thread-MPI threads to start, default=0',
                        type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid",
                        help='List of GPU device id-s to use, default=\"\" ', default="None")
    return parser


if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()

    sys_min = gmx.GmxSys(name=args.name, coor_file=args.f, top_file=args.p)
    sys_min.nt = args.nt
    sys_min.ntmpi = args.ntmpi
    if args.gpuid != "None":
        sys_min.gpu_id = args.gpuid

    sys_min.em_2_steps(out_folder=args.o, name=args.name,
                       no_constr_nsteps=args.min_steps, constr_nsteps=args.min_steps,
                       posres="", create_box_flag=args.box)
    sys_min.convert_trj(traj=False)

    print("\n\nMinimisation was sucessfull \n\tMinimzed directory :\t" + args.o)
    sys_min.display()
