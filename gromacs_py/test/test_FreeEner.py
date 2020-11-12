#!/usr/bin/env python3
# coding: utf-8

"""
Tests for FreeEner class
"""

import pytest
import numpy as np

from gromacs_py.free_ener import FreeEner
from pdb_manip_py import pdb_manip

from .datafiles import PDB_1D30

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_free_bind(tmp_path):
    dna_free = FreeEner('DAP', f'{tmp_path}/1D30_solv')
    dna_free.prepare_complex_pdb(
        PDB_1D30, 'NC(=N)c1ccc(cc1)c2[nH]c3cc(ccc3c2)C(N)=N')

    start_coor = pdb_manip.Coor(dna_free.gmxsys.coor_file)
    assert start_coor.num == 20973
    dna_free.equilibrate_complex(em_steps=100, HA_time=0.001,
                                 CA_time=0.002, CA_LOW_time=0.02,
                                 dt=0.002, dt_HA=0.001, temp=300,
                                 receptor_grp='DNA',
                                 short_steps=500)
    dna_free.compute_add_intermol_from_traj(ref_coor=None, rec_group='DNA')
    dg_rest_water = dna_free.get_water_restr()
    print(dg_rest_water)
    assert ((5.5 < dg_rest_water) and (dg_rest_water < 9.5))

    dna_free.plot_intermol_restr()

    delta_restr = 0.5
    restr_list = np.arange(0, 1 + delta_restr, delta_restr)
    delta_coul = 0.5
    coul_list = np.arange(0, 1 + delta_coul, delta_coul)
    vdw_coul = 0.5
    vdw_list = np.arange(0, 1 + vdw_coul, vdw_coul)
    dna_free.run(lambda_restr_list=restr_list,
                 lambda_coul_list=coul_list,
                 lambda_vdw_list=vdw_list,
                 em_steps=50, nvt_time=0.05, npt_time=0.05, prod_time=0.25,
                 temp_groups='System')
    dna_free.extend_lambda_prod(prod_time=0.5)
    dg_bind, dg_std_bind = dna_free.get_free_ener()
    print(dg_bind, dg_std_bind)
    assert ((0 < dg_bind) and (dg_bind < 200))

    # Need to think about using alchemlyb as a dependance or not
    dna_free.plot_convergence(dt=0.25, graph_out=f'{tmp_path}/alchem.png')
    dna_free.compute_convergence_gbar(dt=0.25)
    dna_free.plot_convergence_graph(graph_out=f'{tmp_path}/bar.png')


@pytest.mark.parametrize("smile, dg", [
    ("c1ccccc1", -1.481),
    ("c1ccccc1O", -0.413),
    ("c1cc(O)ccc1O", -0.826),
])
def test_is_palindrome(smile, dg):
    assert FreeEner.symmetry_correction(smile) == pytest.approx(dg, rel=1e-3)


def test_solv_water_free(tmp_path):
    mol_free = FreeEner('DEN', f'{tmp_path}/indene_water_solv')
    mol_free.water_box_from_SMILE('C1C=Cc2ccccc12')

    start_coor = pdb_manip.Coor(mol_free.gmxsys.coor_file)
    assert start_coor.num == 1625
    # To avoid domain decomposition errors:
    # Use only one proc
    mol_free.gmxsys.nt = 1
    mol_free.equilibrate_solvent_box(em_steps=100, dt=0.002, prod_time=0.02,
                                     short_steps=500)
    delta_coul = 0.5
    coul_list = np.arange(0, 1 + delta_coul, delta_coul)
    vdw_coul = 0.5
    vdw_list = np.arange(0, 1 + vdw_coul, vdw_coul)
    mol_free.run(lambda_coul_list=coul_list,
                 lambda_vdw_list=vdw_list,
                 em_steps=100, nvt_time=0.1, npt_time=0.1, prod_time=0.1,
                 temp_groups='System')
    dg_solv, dg_std_solv = mol_free.get_free_ener()
    print(dg_solv, dg_std_solv)
    assert ((10 < dg_solv) and (dg_solv < 30))


def test_solv_oct(tmp_path):
    mol_free = FreeEner('DEN', f'{tmp_path}/indene_oct_solv')
    mol_free.octanol_box_from_SMILE('C1C=Cc2ccccc12')

    start_coor = pdb_manip.Coor(mol_free.gmxsys.coor_file)
    assert start_coor.num == 1502
    mol_free.gmxsys.nt = 1
    mol_free.equilibrate_solvent_box(em_steps=100, dt=0.002, prod_time=0.02,
                                     short_steps=500)
