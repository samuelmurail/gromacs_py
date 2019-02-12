#!/usr/bin/env python3

""" Minimize a pdb file structure
"""

import argparse
import shutil
import gromacs.gmx5 as gmx

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description='Minimize a cyclic peptide structure in 2 steps, \
        the first step without bonds constraints and the second step with bonds constraints')
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name",
                        help='Output file name', type=str, required=True)
    parser.add_argument('-dir', action="store", dest="out_dir",
                        help='Output directory for intermediate files', type=str, default="./tmp_em")
    parser.add_argument('-m_steps', action="store", dest="min_steps",
                        help='Minimisation nsteps, default=1000', type=int, default=1000)
    parser.add_argument('-keep', action="store_true", dest="keep_flag",
                        help='Flag to keep temporary files (without flag output \
                        directory will be delete')
    parser.add_argument('-cyclic', action="store_true", dest="cyclic_flag",
                        help='Flag to indicate if the peptide/protein is cyclic')
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

    vsite = "none"
    peptide = gmx.GmxSys(name=args.name, coor_file=args.f)
    peptide.nt = args.nt
    peptide.ntmpi = args.ntmpi
    if args.gpuid != "None":
        peptide.gpu_id = args.gpuid

    if args.cyclic_flag:
        peptide.cyclic_peptide_top(out_folder=args.out_dir + '/top')
    else:
        peptide.add_top(out_folder=args.out_dir + '/top',
                        pdb2gmx_option_dict={'vsite': 'no', 'ignh': 'yes', 'ter': 'no'})

    peptide.em(out_folder=args.out_dir + '/em', name="min_" + args.name,
               nsteps=args.min_steps, posres="", nstxout=1000, create_box_flag=True,
               constraints="none")

    peptide.convert_trj(traj=False)

    # Get the minimised structure:
    shutil.copyfile(peptide.coor_file, args.name + '.pdb')

    # Keep or not the intermediate files:
    if not args.keep_flag:
        shutil.rmtree(args.out_dir, ignore_errors=True)
    else:
        peptide.display()
