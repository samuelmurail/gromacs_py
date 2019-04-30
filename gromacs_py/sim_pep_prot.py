#!/usr/bin/env python3

# Minimize a pdb file, return a pdb and a topologie

__author__ = "Samuel Murail"

# Needed because relative imports ..tools don't work
# Need to define package to gromacs_py to import ..tools
# Otherwise package will be gromacs and won't know gromacs_py.tools
#__package__ = 'gromacs_py.gromacs'

import gromacs.gmx5 as gmx
import argparse


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(description="Create a peptide strucure, insert it around a protein and do the minimisation, equilibration and production.")
    # Peptide args:
    parser.add_argument('-seq', action="store", dest="seq", help='Peptide sequence', type=str, required=True)
    parser.add_argument('-npep', action="store", dest="num_mol", help='Number of molecule to insert', type=int, required=True)
    parser.add_argument('-Pep_time', action="store", dest="pep_time", help='Peptide vacuum equilibration with no constraint time(ns), default = 1ns', type=float, default=1)
    # Input sys args:
    parser.add_argument('-fsys', action="store", dest="f_sys", help='Input PDB file of the system', required=True)
    parser.add_argument('-psys', action="store", dest="p_sys", help='Topologie in gromacs format .top of the system', required=True)
    parser.add_argument('-ssys', action="store", dest="s_sys", help='Input tpr file of the system', required=True)
    # General args:
    parser.add_argument('-o', action="store", dest="o", help='Output Directory', required=True)
    parser.add_argument('-n', action="store", dest="name", help='Output file name', required=True)
    parser.add_argument('-dt', action="store", dest="dt", help='Equilibration dt, default=0.005 (5 fs)', type=float, default=0.005)
    parser.add_argument('-em_steps', action="store", dest="em_steps", help='Minimisation steps, default = 5000', type=int, default=5000)
    # Equilibration and production args:
    parser.add_argument('-dt_HA', action="store", dest="dt_HA", help='Equilibration dt, default=0.005 (5 fs)', type=float, default=0.002)
    parser.add_argument('-HA_time', action="store", dest="HA_time", help='Equilibration with HA constraint time(ns), default = 0.25ns', type=float, default=0.25)
    parser.add_argument('-CA_time', action="store", dest="CA_time", help='Equilibration with HA constraint time(ns), default = 1ns', type=float, default=1)
    parser.add_argument('-CA_LOW_time', action="store", dest="CA_LOW_time", help='Equilibration with HA constraint time(ns), default = 5ns', type=float, default=5)
    parser.add_argument('-maxwarn', action="store", dest="maxwarn",
                        help='Total number of warnings allowed for the equilibration, default=0', type=int,
                        default=0)
    parser.add_argument('-PROD_time', action="store", dest="Prod_time", help='Production time(ns), default = 100ns', type=float, default=100)
    # mdrun args:
    parser.add_argument('-nt', action="store", dest="nt", help='Total number of threads to start, default=0', type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi", help='Number of thread-MPI threads to start, default=0', type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid", help='List of GPU device id-s to use, default=\"\" ', default="None")

    return(parser)


if __name__ == "__main__":

    parser = parser_input()
    args = parser.parse_args()


    # General args:
    dt          = args.dt
    out_folder  = args.o
    sys_name    = args.name
    maxwarn     = args.maxwarn

    # Peptide args:
    sequence = args.seq
    em_nsteps = args.em_steps
    pep_step = 1000 * args.pep_time / dt

    # Create peptide:
    peptide = gmx.GmxSys(name='pep_' + sequence)
    peptide.nt = args.nt
    peptide.ntmpi = args.ntmpi
    if args.gpuid != "None":
        peptide.gpu_id = args.gpuid

    peptide.create_peptide(sequence=sequence, out_folder=out_folder + "/" + sequence,
                           em_nsteps=em_nsteps, equi_nsteps=pep_step, posre_post="_pep")
    peptide.display()

    # Starting system args:
    coor_sys = args.f_sys
    top_sys = args.p_sys
    pep_num = args.num_mol

    # Insert peptide:
    sys_pep_prot = gmx.GmxSys(name=sys_name, coor_file=coor_sys, top_file=top_sys, tpr=args.s_sys)
    sys_pep_prot.nt = args.nt
    sys_pep_prot.ntmpi = args.ntmpi
    if args.gpuid != "None":
        sys_pep_prot.gpu_id = args.gpuid

    sys_pep_prot.insert_mol_sys(mol_gromacs=peptide, mol_num=pep_num,
                                new_name=sys_name + "_" + sequence,
                                out_folder=out_folder + "/top_prot_" + sequence + "/")

    # Equilibration and production args
    dt_HA = args.dt_HA
    HA_step = 1000 * args.HA_time / dt_HA
    CA_step = 1000 * args.CA_time / dt
    CA_LOW_step = 1000 * args.CA_LOW_time / dt
    PROD_step = 1000 * args.Prod_time / dt

    sys_pep_prot.em_equi_three_step_iter_error(out_folder=out_folder + "/em_equi_prot_" + sequence + "/", name=sys_name, nsteps_HA=HA_step, nsteps_CA=CA_step, nsteps_CA_LOW=CA_LOW_step, maxwarn=maxwarn)

    sys_pep_prot.production(out_folder=out_folder + "/prod_prot_" + sequence + "/",
                            name=sys_name, nsteps=PROD_step, dt=dt, maxwarn=maxwarn)

    sys_pep_prot.display()
