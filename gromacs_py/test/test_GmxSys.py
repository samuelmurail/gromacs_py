#!/usr/bin/env python3
# coding: utf-8

"""
Tests for GmxSys class
"""

import pytest
import os

from gromacs_py import gmx
import gromacs_py.tools.ambertools as ambertools

from pdb_manip_py import pdb_manip

from .datafiles import PDB_1Y0M, PDB_1RXZ, PDB_1JD4

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_zinc_finger(tmp_path):

    gmx.show_log()
    ##################################
    # ##   Create the topologie:   ###
    ##################################
    prot = gmx.GmxSys(name='1jd4', coor_file=PDB_1JD4)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_DIAP1'),
                     vsite='hydrogens', ignore_ZN=False)

    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 3292

    prot.switch_ion_octa_dummy()
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 3304

    prot.solvate_add_ions(out_folder=os.path.join(tmp_path, 'top_sys'),
                          maxwarn=3)

    prot.em(out_folder=os.path.join(tmp_path, 'em_DIAP1'),
            nsteps=10,
            constraints='none')
    em_coor = pdb_manip.Coor(prot.coor_file)
    assert em_coor.num == 56366

    ener_pd = prot.get_ener(['Potential'])
    assert ener_pd['Potential'].values[-1] < -155000


def test_set_protonate_charmm(tmp_path):

    gmx.show_log()
    ##########################################
    # ##   Create the topologie at PH 4.0  ###
    ##########################################

    res_prot_dict = {'GLUP': {'res_num': [9, 24, 31, 57]},
                     'ASPP': {'res_num': [3, 70]},
                     'HSP': {'res_num': [28, 35]},
                     'HSD': {'res_num': [175]},
                     'HSE': {'res_num': [214]}}

    prot = gmx.GmxSys(name='1RXZ_ph4', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph4'), ph=7.0,
                     res_prot_dict=res_prot_dict)
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4081

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -4


def test_set_protonate_amber(tmp_path):

    gmx.show_log()
    ##########################################
    # ##   Create the topologie at PH 4.0  ###
    ##########################################

    res_prot_dict = {'GLH': {'res_num': [9, 24, 31, 57]},
                     'ASH': {'res_num': [3, 70]},
                     'HIP': {'res_num': [28, 35]},
                     'HID': {'res_num': [175]},
                     'HIE': {'res_num': [214]}}

    prot = gmx.GmxSys(name='1RXZ_ph4', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph4'), ph=7.0,
                     res_prot_dict=res_prot_dict,
                     ff="amber99sb-ildn")
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4081

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -4


def test_protonate_ph4_7_12_charmm(tmp_path):

    gmx.show_log()
    ##########################################
    # ##   Create the topologie at PH 4.0  ###
    ##########################################
    prot = gmx.GmxSys(name='1RXZ_ph4', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph4'), ph=4.0)
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4108

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == 23

    ##########################################
    # ##   Create the topologie at PH 7.0  ###
    ##########################################
    prot = gmx.GmxSys(name='1RXZ_ph7', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph7'), ph=7.0)
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4073

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -12

    ##########################################
    # ##   Create the topologie at PH 10.0 ###
    ##########################################

    # Remark: Charmm doesn't allow lysine unprotonated
    prot = gmx.GmxSys(name='1RXZ_ph10', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph10'),
                     ph=10.0)
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4073

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -12


def test_protonate_ph4_7_amber(tmp_path):

    gmx.show_log()
    ##########################################
    # ##   Create the topologie at PH 4.0  ###
    ##########################################
    prot = gmx.GmxSys(name='1RXZ_ph4', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph4'),
                     ph=4.0, ff="amber99sb-ildn")
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4108

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == 23

    ##########################################
    # ##   Create the topologie at PH 7.0  ###
    ##########################################
    prot = gmx.GmxSys(name='1RXZ_ph7', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph7'),
                     ph=7.0, ff="amber99sb-ildn")
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4073

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -12

    ##########################################
    # ##   Create the topologie at PH 10.0 ###
    ##########################################
    # Remark: Amber allow lysine unprotonated
    prot = gmx.GmxSys(name='1RXZ_ph10', coor_file=PDB_1RXZ)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_1RXZ_ph10'),
                     ph=10.0, ff="amber99sb-ildn")
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 4071

    top_1RXZ = gmx.TopSys(prot.top_file)
    assert top_1RXZ.charge() == -14


def test_insert_ethanol(tmp_path):

    gmx.show_log()
    ##################################
    # ##   Create the topologie:   ###
    ##################################
    prot = gmx.GmxSys(name='1y0m', coor_file=PDB_1Y0M)
    prot.prepare_top(out_folder=os.path.join(tmp_path, 'top_SH3'))
    top_coor = pdb_manip.Coor(prot.coor_file)
    assert top_coor.num == 996

    prot.solvate_add_ions(out_folder=os.path.join(tmp_path, 'top_sys'))

    prot.em(out_folder=os.path.join(tmp_path, 'em_SH3'),
            nsteps=10,
            constraints='none')
    em_coor = pdb_manip.Coor(prot.coor_file)
    assert em_coor.num == 15383

    ener_pd = prot.get_ener(['Potential'])
    assert ener_pd['Potential'].values[-1] < -155000

    ###################################
    # ###   Create Ethanol system   ###
    ###################################
    smile = 'CCO'
    eth_pdb = os.path.join(tmp_path, 'ethanol.pdb')
    charge = ambertools.smile_to_pdb(smile, eth_pdb, 'ETH')

    assert charge == 0

    eth_sys = gmx.GmxSys(name='ETH', coor_file=eth_pdb)
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


@pytest.mark.skipif(float(gmx.gmxsys.gmx_version) >= 2019,
                    reason="Gromacs verions >= 19 have issues with ACE")
def test_insert_peptide_vsite(tmp_path):

    gmx.show_log()
    ##################################
    # ##   Create the topologie:   ###
    ##################################
    prot = gmx.GmxSys(name='1y0m', coor_file=PDB_1Y0M)
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
    # ##    Create a DD peptide     ###
    ###################################
    pep = gmx.GmxSys(name='DD')
    pep.create_peptide(sequence='DD',
                       out_folder=os.path.join(tmp_path, 'top_DD'),
                       em_nsteps=10, equi_nsteps=0,
                       vsite='hydrogens')

    ####################################################
    # # Insert 4 copy of ethanol in  the SH3 system: ###
    ####################################################

    prot.insert_mol_sys(mol_gromacs=pep, mol_num=10,
                        new_name='SH3_DD',
                        out_folder=os.path.join(tmp_path, 'top_DD_SH3'))
    ################################
    # ##   Minimize the system   ###
    ################################
    prot.em_2_steps(out_folder=os.path.join(tmp_path, 'em_DD_SH3'),
                    no_constr_nsteps=10, constr_nsteps=10)
    prot_lig_coor = pdb_manip.Coor(prot.coor_file)
    assert 15300 < prot_lig_coor.num < 15518

    prot.display_history()
