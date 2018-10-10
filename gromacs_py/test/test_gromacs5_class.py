#!/usr/bin/env python3

##################################
#########   TEST PART   ##########
##################################

__author__ = "Samuel Murail"

import sys
import os

# necessary to load modules in root dir 
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import gromacs.gmx5 as gmx
import tools.pdb_manip as pdb
import tools.osCommand as osCommand



if __name__ == '__main__':


    dir_out = "./output/"
   
    test_sys = gmx.gmx_sys(name = "test", coor_file = "input/1AWR-055_bestene1-mc.pdb")
    test_sys.prepare_top(out_folder = dir_out+"/tmp/prot/top", vsite = "hydrogens")
    test_sys.create_box()
    test_sys.em_2_steps(out_folder = dir_out+"/tmp/prot/em", no_constr_nsteps = 500, constr_nsteps = 500)

    test_sys.solvate_add_ions(out_folder = dir_out+"/tmp/sys/top")
    test_sys.em_2_steps(out_folder = dir_out+"/tmp/sys/em", no_constr_nsteps = 500, constr_nsteps = 500)
    test_sys.equi_three_step(out_folder = dir_out+"/tmp/sys/equi", nsteps_HA = 1000, nsteps_CA = 1000, nsteps_CA_LOW = 1000)
    test_sys.production(out_folder = dir_out+"/tmp/sys/prod", nsteps = 1000)
    #test_sys.extend_equi_prod(out_folder = dir_out+"/tmp/sys/prod", nsteps = 10000)
    test_sys.display()

    #del test_sys
    test_sys_2 = gmx.gmx_sys(name = "test", coor_file = "output/tmp/sys/prod/prod_test.gro")
    test_sys_2.display()


    peptide = gmx.gmx_sys(name = "peptide_AP")
    peptide.create_peptide(sequence = "AP", out_folder = "output/tmp/AP/top/")

    test_sys.insert_mol_sys(mol_gromacs = peptide, mol_num = 10, new_name = "test_AP", out_folder = "output/tmp/sys_AP/top/")
    test_sys.em_2_steps(out_folder = "output/tmp/sys_AP/em/", no_constr_nsteps = 500, constr_nsteps = 500)
    test_sys.convert_trj(traj = False)    
    test_sys.display()

    #dir_out = "./output/"
#
    #GPH = gmx.create_peptide(sequence = "APH", out_folder = dir_out+"APH", em_nsteps = 1000, equi_nsteps = 1000)
    #print(GPH)
#
#
    ## Topologie creation:
    #prot = gmx.prepare_top(pdb_in = "input/1AWR-055_bestene1-mc.pdb",
    #    out_folder = dir_out+"1AWR-055/00_param", sys_name = "1AWR_055")
#
    ## Minimize the protein :
    #prot_min = gmx.em_2_steps(pdb_in = prot['f'],  top_in = prot['p'], out_folder = dir_out+"1AWR-055/00_param/em_prot",
    #    sys_name = "mini_1AWR_055", no_constr_nsteps = 500, constr_nsteps = 500, posres = "-DPOSRES_HA_LOW", create_box_flag = True)
#
    ## Solvate protein :
    #CYPA = gmx.solvate_add_ions(pdb_in = prot_min['f'], top_in = prot_min['p'], out_folder = dir_out+"1AWR-055/00_param/",
    #    sys_name = "1AWR_055")
#
    ## Equilibrate and production:
    #CYPA_mini = gmx.em_2_steps(pdb_in = CYPA['f'], top_in = CYPA['p'],
    #    out_folder = dir_out+"1AWR-055/01_mini/", sys_name = "1AWR_055",
    #    no_constr_nsteps = 500, constr_nsteps = 500)
    #CYPA_equi = gmx.equi_three_step(pdb_in = CYPA_mini['f'], top_in = CYPA_mini['p'],
    #    out_folder = dir_out+"1AWR-055/02_equi/", sys_name = "1AWR_055",
    #    nsteps_HA = 1000, nsteps_CA = 1000, nsteps_CA_LOW = 1000)
    #CYPA_prod = gmx.production(pdb_in = CYPA_equi['f'], top_in = CYPA_equi['p'],
    #    out_folder= dir_out+"1AWR-055/03_prod/", sys_name = "1AWR_055", nsteps = 100)
#
    #print(CYPA_prod)
#
    #CYPA_prod_compact = gmx.convert_trj(f = CYPA_prod['f'], o = CYPA_prod['f'][:-4]+"_compact.pdb", s = CYPA_prod['f'][:-4]+".tpr",
    #    ur = "compact", pbc = "mol", select = "System", check_file_out = False)
#
    #pep_prot_sys = gmx.insert_mol_sys(sys_pdb = CYPA_prod_compact , sys_top = CYPA_prod['p'],
    #    mol_pdb = GPH['f'], mol_top = GPH['p'], mol_num = 20, sys_name = "CYPA_20_GPH",
    #    out_folder = dir_out+"CYPA_20_GPH")
    ##mol_length = 3,
    #Cypa_GPH_min = gmx.em_2_steps(pdb_in = pep_prot_sys['f'], top_in = pep_prot_sys['p'], out_folder = dir_out+"CYPA_20_GPH/01_mini_no_constraints/",
    #    sys_name = "CYPA_20_GPH", create_box_flag = False, no_constr_nsteps = 500, constr_nsteps = 500)
#
#
    #CYPA_equi = gmx.equi_three_step(pdb_in = Cypa_GPH_min['f'], top_in = Cypa_GPH_min['p'],
    #    out_folder = dir_out+"CYPA_20_GPH/02_mini_equi_3_steps/", sys_name = "CYPA_20_GPH",
    #    nsteps_HA = 1000, nsteps_CA = 1000, nsteps_CA_LOW = 1000)
