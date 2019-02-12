#!/usr/bin/env python3

__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Simulation production")
    parser.add_argument('-f', action="store", dest="f", help='Input PDB file', type=str, required=True)
    parser.add_argument('-p', action="store", dest="p", help='Topologie in gromacs format .top', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o", help='Output Directory', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name", help='Output file name', type=str, required=True)
    parser.add_argument('-time', action="store", dest="time", help='Production time, default=10', type=float, default=10)
    parser.add_argument('-dt', action="store", dest="dt", help='Equilibration dt, default=0.005 (5 fs)', type=float, default=0.005)
    parser.add_argument('-nt', action="store", dest="nt", help='Total number of threads to start, default=0', type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi", help='Number of thread-MPI threads to start, default=0', type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid", help='List of GPU device id-s to use, default=\"\" ', default="None")
    return(parser)


if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()

    dt = args.dt
    nsteps = 1000 * args.time / dt

    sys_prod = gmx.GmxSys(name=args.name, coor_file=args.f, top_file=args.p)
    sys_prod.nt = args.nt
    sys_prod.ntmpi = args.ntmpi
    if args.gpuid != "None":
        sys_prod.gpu_id = args.gpuid

    sys_prod.production(out_folder=args.o, name=args.name, nsteps=nsteps, dt=args.dt)

    print("\n\nProductuion was sucessfull \n\tProduction directory :\t" + args.o)

    sys_prod.display()
