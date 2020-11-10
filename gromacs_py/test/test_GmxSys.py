#!/usr/bin/env python3
# coding: utf-8

"""
Tests for GmxSys class
"""

import pytest
import numpy as np
import os

import gromacs_py.gmx
from gromacs_py.gmx import GmxSys
import gromacs_py.tools.ambertools as ambertools

from pdb_manip_py import pdb_manip

from .datafiles import PDB_1Y0M

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_insert_mol(tmp_path):

    gromacs_py.gmx.show_log()
    ##################################
    # ##   Create the topologie:   ###
    ##################################
    prot = GmxSys(name='1y0m', coor_file=PDB_1Y0M)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_SH3'),
                     vsite='hydrogens')
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 1064

    prot.solvate_add_ions(out_folder=os.path.join(tmp_path, 'top_sys'))

    prot.em(out_folder=os.path.join(tmp_path, 'em_SH3'),
            nsteps=10,
            constraints='none')
    em_coor = pdb_manip.Coor(prot.coor_file)
    assert em_coor.num == 15337

    ener_pd = prot.get_ener(['Potential'])
    assert ener_pd['Potential'].values[-1] < -155000

    ###################################
    # ###   Create Ethanol system   ###
    ###################################
    smile = 'CCO'
    eth_pdb = os.path.join(tmp_path, 'ethanol.pdb')
    charge = ambertools.smile_to_pdb(smile, eth_pdb, 'ETH')

    assert charge == 0

    eth_sys = GmxSys(name='ETH', coor_file=eth_pdb)
    eth_sys.prepare_top_ligand(
        out_folder=os.path.join(tmp_path, 'eth_top'),
        ff='amber99sb-ildn', include_mol={'ETH': smile})

    eth_sys.create_box(name=None, dist=0.5)

    ####################################################
    # # Insert 4 copy of ethanol in  the SH3 system: ###
    ####################################################

    prot.insert_mol_sys(mol_gromacs=eth_sys, mol_num=10,
                        new_name='SH3_ETH',
                        out_folder=os.path.join(tmp_path, 'top_ETH_SH3'))
    ################################
    # ##   Minimize the system   ###
    ################################
    prot.em_2_steps(out_folder=os.path.join(tmp_path, 'em_ETH_SH3'),
                    no_constr_nsteps=10, constr_nsteps=10)
    prot_lig_coor = pdb_manip.Coor(prot.coor_file)
    assert 15300 < prot_lig_coor.num < 15518

    prot.display_history()

    """
    >>> ################################
    >>> ####   Minimize the system   ###
    >>> ################################
    >>> prot.em_2_steps(out_folder=os.path.join(tmp_path, 'top_D_SH3'), \
no_constr_nsteps=10, constr_nsteps=10)
    - Create the tpr file Init_em_1y0m.tpr
    gmx grompp -f Init_em_1y0m.mdp -c SH3_D_neutral.gro -r SH3_D_neutral.gro \
-p SH3_D_neutral.top -po out_Init_em_1y0m.mdp -o Init_em_1y0m.tpr -maxwarn 1
    - Launch the simulation Init_em_1y0m.tpr
    gmx mdrun -s Init_em_1y0m.tpr -deffnm Init_em_1y0m -nt 0 -ntmpi 0 \
-nsteps -2 -nocopyright
    - Create the tpr file 1y0m.tpr
    gmx grompp -f 1y0m.mdp -c Init_em_1y0m.gro -r Init_em_1y0m.gro -p \
SH3_D_neutral.top -po out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
    - Launch the simulation 1y0m.tpr
    gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 -nocopyright
    >>> ##################################
    >>> ####    Show system history    ###
    >>> ##################################
    >>> prot.display_history() #doctest: +ELLIPSIS
    State -3:
    <BLANKLINE>
    name         : 1y0m
    sim_name     : genion_1y0m_water_ion
    coor_file    : .../top_sys/1y0m_water_ion.gro
    top_file     : .../top_sys/1y0m_water_ion.top
    tpr          : .../top_sys/genion_1y0m_water_ion.tpr
    mdp          : ...template/mini.mdp
    nt           : 0
    ntmpi        : 0
    sys_history  : 0
    <BLANKLINE>
    State -2:
    <BLANKLINE>
    name         : 1y0m
    sim_name     : genion_SH3_D_neutral
    coor_file    : .../top_D_SH3/SH3_D_neutral.gro
    top_file     : .../top_D_SH3/SH3_D_neutral.top
    tpr          : .../top_D_SH3/genion_SH3_D_neutral.tpr
    mdp          : ...template/mini.mdp
    xtc          : .../em_SH3/1y0m.trr
    edr          : .../em_SH3/1y0m.edr
    log          : .../em_SH3/1y0m.log
    nt           : 0
    ntmpi        : 0
    sys_history  : 0
    <BLANKLINE>
    State -1:
    <BLANKLINE>
    name         : 1y0m
    sim_name     : Init_em_1y0m
    coor_file    : .../top_D_SH3/Init_em_1y0m.gro
    top_file     : .../top_D_SH3/SH3_D_neutral.top
    tpr          : .../top_D_SH3/Init_em_1y0m.tpr
    mdp          : .../top_D_SH3/Init_em_1y0m.mdp
    xtc          : .../top_D_SH3/Init_em_1y0m.trr
    edr          : .../top_D_SH3/Init_em_1y0m.edr
    log          : .../top_D_SH3/Init_em_1y0m.log
    nt           : 0
    ntmpi        : 0
    sys_history  : 0
    <BLANKLINE>
    >>> ###################################
    >>> ####   Equilibrate the system   ###
    >>> ###################################
    >>> equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME, \
"template/equi_vsites.mdp")
    >>> mdp_options = {'nsteps': 100, 'define': '-DPOSRES', 'dt': 0.001}
    >>> prot.run_md_sim(out_folder=os.path.join(tmp_path, 'equi_HA_D_SH3'), \
name="equi_HA_D_SH3", mdp_template=equi_template_mdp,\
                        mdp_options=mdp_options)
    - Create the tpr file equi_HA_D_SH3.tpr
    gmx grompp -f equi_HA_D_SH3.mdp -c ../top_D_SH3/1y0m.gro -r \
../top_D_SH3/1y0m.gro -p ../top_D_SH3/SH3_D_neutral.top -po \
out_equi_HA_D_SH3.mdp -o equi_HA_D_SH3.tpr -maxwarn 0
    - Launch the simulation equi_HA_D_SH3.tpr
    gmx mdrun -s equi_HA_D_SH3.tpr -deffnm equi_HA_D_SH3 -nt 0 -ntmpi 0 \
-nsteps -2 -nocopyright
    >>> prot.get_simulation_time() #doctest: +ELLIPSIS
    - Get simulation time from : .../equi_HA_D_SH3/equi_HA_D_SH3.cpt
    gmx check -f .../equi_HA_D_SH3/equi_HA_D_SH3.cpt
    0.1
    >>> prot.convert_trj(traj=False) #doctest: +ELLIPSIS
    - Convert trj/coor
    gmx trjconv -f .../equi_HA_D_SH3/equi_HA_D_SH3.gro -o \
.../equi_HA_D_SH3/equi_HA_D_SH3_compact.pdb -s \
.../equi_HA_D_SH3/equi_HA_D_SH3.tpr -ur compact -pbc mol
    >>> prot.display() #doctest: +ELLIPSIS
    name         : 1y0m
    sim_name     : equi_HA_D_SH3
    coor_file    : .../equi_HA_D_SH3/equi_HA_D_SH3_compact.pdb
    top_file     : .../top_D_SH3/SH3_D_neutral.top
    tpr          : .../equi_HA_D_SH3/equi_HA_D_SH3.tpr
    mdp          : .../equi_HA_D_SH3/equi_HA_D_SH3.mdp
    xtc          : .../equi_HA_D_SH3/equi_HA_D_SH3.xtc
    edr          : .../equi_HA_D_SH3/equi_HA_D_SH3.edr
    log          : .../equi_HA_D_SH3/equi_HA_D_SH3.log
    nt           : 0
    ntmpi        : 0
    sys_history  : 4
    >>> #########################################
    >>> ### Extract Potential Energy and Temp ###
    >>> #########################################
    >>> ener_pd = prot.get_ener(['Potential', 'Temp'])  #doctest: +ELLIPSIS
    - Extract energy
    gmx energy -f .../equi_HA_D_SH3/equi_HA_D_SH3.edr -o tmp_edr.xvg
    >>> ener_pd['Potential'].mean() #doctest: +ELLIPSIS
    -2...
    >>> rmsd_pd = prot.get_rmsd(['C-alpha', 'Protein'])  #doctest: +ELLIPSIS
    - Extract RMSD
    - Create the ndx file ...equi_HA_D_SH3.ndx
    gmx make_ndx -f ...equi_HA_D_SH3_compact.pdb -o ...equi_HA_D_SH3.ndx
    gmx rms -s ...equi_HA_D_SH3.tpr -f ...equi_HA_D_SH3.xtc -n ...\
equi_HA_D_SH3.ndx -o tmp_rmsd.xvg -fit rot+trans -ng 1 -pbc no
    >>> rmsd_pd #doctest: +ELLIPSIS
       time ...Protein
    0   0.0...
    >>> rmsf_pd = prot.get_rmsf(['Protein'], res="yes")  #doctest: +ELLIPSIS
    - Extract RMSF
    gmx rmsf -s ...equi_HA_D_SH3.tpr -f ...equi_HA_D_SH3.xtc -n ...\
equi_HA_D_SH3.ndx -o tmp_rmsf.xvg -fit no -res yes
    >>> rmsf_pd #doctest: +ELLIPSIS
        Residue    RMSF
    0       791  ...
    """