#!/usr/bin/env python3
# coding: utf-8

"""
Tests for GmxSys class
"""

import os

from gromacs_py import gmx

from pdb_manip_py import pdb_manip

from .datafiles import PDB_5VAV

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_insert_peptide_vsite(tmp_path):
    # Gromacs verions >= 2021 support cyclic peptides

    gmx.show_log()

    cyclic_pep = gmx.GmxSys(
        name='5vav',
        coor_file=PDB_5VAV)

    top_coor = pdb_manip.Multi_Coor(cyclic_pep.coor_file)
    assert top_coor.coor_list[0].num == 209

    cyclic_pep.cyclic_peptide_top(
        out_folder=os.path.join(tmp_path, 'cyclic/top'))

    top_coor = pdb_manip.Coor(cyclic_pep.coor_file)
    assert top_coor.num == 209

    cyclic_top = gmx.TopSys(cyclic_pep.top_file)

    assert cyclic_top.prot_res_num() == 14

    cyclic_top.display()

    cyclic_pep.em(
        out_folder=os.path.join(tmp_path, 'cyclic/em/'),
        nsteps=10,
        create_box_flag=True)

    cyclic_amber_pep = gmx.GmxSys(
        name='5vav_amber',
        coor_file=PDB_5VAV)

    top_coor = pdb_manip.Coor(cyclic_pep.coor_file)
    assert top_coor.num == 209

    cyclic_amber_pep.cyclic_peptide_top(
        out_folder=os.path.join(tmp_path, 'cyclic/top'),
        ff='amber99sb-ildn')

    top_coor = pdb_manip.Coor(cyclic_pep.coor_file)
    assert top_coor.num == 209

    cyclic_amber_top = gmx.TopSys(cyclic_amber_pep.top_file)

    print(cyclic_amber_top.charge())

    cyclic_amber_pep.em(
        out_folder=os.path.join(tmp_path, 'cyclic/em/'),
        nsteps=10,
        create_box_flag=True)
