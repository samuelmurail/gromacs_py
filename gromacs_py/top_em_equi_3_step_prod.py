#!/usr/bin/env python3


__author__ = "Samuel Murail"


import gromacs.gmx5 as gmx
import argparse


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(
        description="Equilibrate in 3 steps a system (coor+top), (i) first equilibration with heavy atoms position restraints, (ii)"
        " second equilibration with alpha carbon position restraints and (iii) finaly equilibration with weak alpha carbon position restraints")
    parser.add_argument('-f', action="store", dest="f", help='Input PDB file', type=str, required=True)
    parser.add_argument('-o', action="store", dest="o", help='Output Directory', type=str, required=True)
    parser.add_argument('-vsite', action="store_true", dest="vsite_flag", help='Use virtual site for hydrogens')
    parser.add_argument('-C', action="store", dest="Conc", help='Ion concentration (mM), default = 0.15 (150mM)', type=float, default=0.15)
    parser.add_argument('-m_steps', action="store", dest="min_steps", help='Minimisation nsteps, default=1000', type=int, default=1000)
    parser.add_argument('-HA_time', action="store", dest="HA_time", help='Equilibration with HA constraint time(ns), default = 0.25ns', type=float, default=0.25)
    parser.add_argument('-CA_time', action="store", dest="CA_time", help='Equilibration with HA constraint time(ns), default = 1ns', type=float, default=1)
    parser.add_argument('-CA_LOW_time', action="store", dest="CA_LOW_time", help='Equilibration with HA constraint time(ns), default = 5ns', type=float, default=5)
    parser.add_argument('-prod_time', action="store", dest="prod_time", help='Production time, default=10', type=float, default=10)
    parser.add_argument('-dt_HA', action="store", dest="dt_HA", help='Equi HA dt, default=0.005 (5 fs)', type=float, default=0.002)
    parser.add_argument('-dt', action="store", dest="dt", help='Equi CA, CA_LOW, dt, default=0.005 (5 fs)', type=float, default=0.005)
    parser.add_argument('-nt', action="store", dest="nt", help='Total number of threads to start, default=0', type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi", help='Number of thread-MPI threads to start, default=0', type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid", help='List of GPU device id-s to use, default=\"\" ', default="None")

    return(parser)


if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()

    if args.vsite_flag:
        vsite = "hydrogens"
    else:
        vsite = "none"

    dt_HA = args.dt_HA
    dt = args.dt
    HA_step = 1000 * args.HA_time / dt_HA
    CA_step = 1000 * args.CA_time / dt
    CA_LOW_step = 1000 * args.CA_LOW_time / dt
    prod_steps = 1000 * args.prod_time / dt
    sys_name = args.f.split("/")[-1][:-4]

    print("EM steps : {}\nEqui HA time : {}ns \nEqui CA time : {}ns \n Equi CA_LOW time :{}ns \n".
          format(args.min_steps, args.HA_time, args.CA_time, args.CA_LOW_time, args.prod_time))

    top_dir = args.o + "/prot_top/"

    md_sys = gmx.GmxSys(name=sys_name, coor_file=args.f)
    md_sys.nt = args.nt
    md_sys.ntmpi = args.ntmpi
    if args.gpuid != "None":
        md_sys.gpu_id = args.gpuid

    # TOPOLOGIE
    md_sys.prepare_top(out_folder=top_dir, vsite=vsite)
    md_sys.create_box(dist=1.0, box_type="dodecahedron", check_file_out=True)
    print("\n\nTopologie creation was sucessfull \n\tTopologie directorie :\t" + top_dir)
    md_sys.display()

    # EM
    em_dir = args.o + "/prot_em/"
    md_sys.em_2_steps(out_folder=em_dir, no_constr_nsteps=args.min_steps, constr_nsteps=args.min_steps,
                      posres="", create_box_flag=False)
    print("\n\nMinimisation was sucessfull \n\tMinimzed directory :\t" + em_dir)
    md_sys.display()

    # SOLVATION
    solv_dir = args.o + "/sys_top/"
    md_sys.solvate_add_ions(out_folder=solv_dir, name=sys_name, ion_C=args.Conc, create_box_flag=False)
    print("\n\nSolvation was sucessfull \n\tTopologie directorie :\t" + solv_dir)
    md_sys.display()

    # EM and  EQUI
    equi_dir = args.o + "/sys_em_equi/"
    md_sys.em_equi_three_step_iter_error(out_folder=equi_dir, name=sys_name, nsteps_HA=HA_step,
                           nsteps_CA=CA_step, nsteps_CA_LOW=CA_LOW_step, dt=dt, dt_HA=dt_HA)
    print("\n\nEquilibration was sucessfull \n\tEquilibration directory :\t" + equi_dir)
    md_sys.display()

    # PROD
    prod_dir = args.o + "/sys_prod/"
    md_sys.production(out_folder=prod_dir, name=sys_name, nsteps=prod_steps, dt=dt)
    print("\n\nProductuion was sucessfull \n\tProduction directory :\t" + prod_dir)
    md_sys.display()
