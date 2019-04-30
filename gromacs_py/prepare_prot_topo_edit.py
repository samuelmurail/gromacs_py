#!/usr/bin/env python3

""" Equilibrate in 3 steps a system (coor+top), (i) first equilibration with heavy \
atoms position restraints, (ii) second equilibration with alpha carbon position \
restraints and (iii) finaly equilibration with weak alpha carbon position restraints
"""

import argparse
import gromacs.gmx5 as gmx

print(gmx.__file__)

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="(i) Create topologie for a protein, (ii) minimize the protein structure in 2 steps,\
                                     (iii) Sovate and add ions, (iv) minimize the system structure, \
                                     (v) first equilibration with heavy atoms position restraints, (vi) second equilibration with \
                                     alpha carbon position restraints and (vii) finaly equilibration with weak alpha carbon \
                                     position restraints")
    # Input
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', required=True)
    parser.add_argument('-p', action="store",dest="p",
                        help='Input topology file', default="None")
    # Output
    parser.add_argument('-o', action="store", dest="o",
                        help='Output Directory', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name",
                        help='Output file name', type=str, required=True)

    # Options
    parser.add_argument('-no_vsite', action="store_true", dest="novsite_flag",
                        help='Use virtual site for hydrogens')
    parser.add_argument('-C', action="store", dest="Conc",
                        help='Ion concentration (mM), default = 0.15 (150mM)',
                        type=float, default=0.15)

    parser.add_argument('-m_steps', action="store", dest="min_steps",
                        help='Minimisation nsteps, default=10000', type=int, default=10000)
    parser.add_argument('-HA_time', action="store", dest="HA_time",
                        help='Equilibration with HA constraint time(ns), default = 2.5 ns',
                        type=float, default=2.5)
    parser.add_argument('-CA_time', action="store", dest="CA_time",
                        help='Equilibration with HA constraint time(ns), default = 5 ns',
                        type=float, default=5)
    parser.add_argument('-CA_LOW_time', action="store", dest="CA_LOW_time",
                        help='Equilibration with HA constraint time(ns), default = 10 ns',
                        type=float, default=10)
    parser.add_argument('-maxwarn', action="store", dest="maxwarn",
                        help='Total number of warnings allowed for the equilibration, default=0', type=int,
                        default=0)
    parser.add_argument('-dt_HA', action="store", dest="dt_HA",
                        help='Equi HA dt, default=0.002 (2 fs)', type=float, default=0.002)
    parser.add_argument('-dt', action="store", dest="dt",
                        help='Equi CA, CA_LOW, dt, default=0.005 (5 fs)', type=float,
                        default=0.005)
    parser.add_argument('-nt', action="store", dest="nt",
                        help='Total number of threads to start, default=0', type=float,
                        default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi",
                        help='Number of thread-MPI threads to start, default=0', type=float,
                        default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid",
                        help='List of GPU device id-s to use, default=\"\" ', default="None")

    return parser


if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()

    if not args.novsite_flag:
        vsite = "hydrogens"
    else:
        vsite = "none"
        if (args.dt_HA > 0.002) or (args.dt > 0.002):
            print("Wrong dt, using dt >0.002 with vsite")
            exit

    sys_name = args.name
    maxwarn = args.maxwarn
    min_steps = args.min_steps
    dt_HA = args.dt_HA
    dt = args.dt
    HA_step = 1000 * args.HA_time / dt_HA
    CA_step = 1000 * args.CA_time / dt
    CA_LOW_step = 1000 * args.CA_LOW_time / dt

    prot_top_folder = args.o + "/top_prot/"
    prot_min_folder = args.o + "/em_prot/"

    sys_top_folder = args.o + "/top_sys/"
    sys_em_equi_folder = args.o + "/em_equi_sys/"

    if args.p != "None":
        prot_sys = gmx.GmxSys(name=sys_name, coor_file=args.f, top_file=args.p)
    else:
        prot_sys = gmx.GmxSys(name=sys_name, coor_file=args.f)
    prot_sys.nt = args.nt
    prot_sys.ntmpi = args.ntmpi
    if args.gpuid != "None":
        prot_sys.gpu_id = args.gpuid

    if args.p == "None":
        prot_sys.prepare_top(out_folder=prot_top_folder, vsite=vsite)

    prot_sys.create_box(dist=1.0, box_type="dodecahedron", check_file_out=True)

    prot_sys.em_2_steps(out_folder=prot_min_folder, name=sys_name,
                        no_constr_nsteps=min_steps, constr_nsteps=min_steps,
                        posres="", create_box_flag=False)

    prot_sys.convert_trj(traj=False)

    prot_sys.solvate_add_ions(out_folder=sys_top_folder, name=sys_name, ion_C=args.Conc)

    prot_sys.em_equi_three_step_iter_error(out_folder=sys_em_equi_folder, name=sys_name, nsteps_HA=HA_step, nsteps_CA=CA_step, nsteps_CA_LOW=CA_LOW_step, maxwarn=maxwarn)

    prot_sys.convert_trj(traj=False)

    print("\n\nSystem preparation was sucessfull \n\tOutput directory :\t" + args.o)

    prot_sys.display()
