#!/usr/bin/env python3

""" Equilibrate in 3 steps a system (coor+top), (i) first equilibration with heavy \
atoms position restraints, (ii) second equilibration with alpha carbon position \
restraints and (iii) finaly equilibration with weak alpha carbon position restraints
"""
import argparse
import gromacs.gmx5 as gmx

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Equilibrate in 3 steps a system (coor+top), (i) \
        first equilibration with heavy atoms position restraints, (ii) second equilibration with \
        alpha carbon position restraints and (iii) finaly equilibration with weak alpha carbon \
        position restraints")
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', type=str, required=True)
    parser.add_argument('-p', action="store", dest="p",
                        help='Topologie in gromacs format .top', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o",
                        help='Output Directory', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name",
                        help='Output file name', type=str, required=True)
    parser.add_argument('-HA_time', action="store", dest="HA_time",
                        help='Equilibration with HA constraint time(ns), default = 0.25ns',
                        type=float, default=0.25)
    parser.add_argument('-CA_time', action="store", dest="CA_time",
                        help='Equilibration with HA constraint time(ns), default = 1ns',
                        type=float, default=1)
    parser.add_argument('-CA_LOW_time', action="store", dest="CA_LOW_time",
                        help='Equilibration with HA constraint time(ns), default = 5ns',
                        type=float, default=5)
    parser.add_argument('-dt_HA', action="store", dest="dt_HA",
                        help='Equi HA dt, default=0.005 (5 fs)', type=float, default=0.002)
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

    dt_HA = args.dt_HA
    dt = args.dt
    maxwarn = args.maxwarn
    HA_step = 1000 * args.HA_time / dt_HA
    CA_step = 1000 * args.CA_time / dt
    CA_LOW_step = 1000 * args.CA_LOW_time / dt

    print("\nEqui HA time :", args.HA_time,
          "ns\nEqui CA time :", args.CA_time,
          "ns\nEqui CA_LOW time :", args.CA_LOW_time, "ns")

    sys_equi = gmx.GmxSys(name=args.name, coor_file=args.f, top_file=args.p)
    sys_equi.nt = args.nt
    sys_equi.ntmpi = args.ntmpi
    if args.gpuid != "None":
        sys_equi.gpu_id = args.gpuid

    sys_equi.equi_three_step(out_folder=args.o, name=args.name, nsteps_HA=HA_step,
                             nsteps_CA=CA_step, nsteps_CA_LOW=CA_LOW_step, dt=dt, dt_HA=dt_HA)

    print("\n\nEquilibration was sucessfull \n\tEquilibration directory :\t" + args.o)

    sys_equi.display()
