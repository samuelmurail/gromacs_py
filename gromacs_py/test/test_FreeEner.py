import pytest
import numpy as np

from gromacs_py.free_ener import FreeEner
from pdb_manip_py import pdb_manip

# Need to add pytest out folder 


@pytest.mark.parametrize("smile, dg", [
    ("c1ccccc1", -1.481),
    ("c1ccccc1O", -0.413),
    ("c1cc(O)ccc1O", -0.826),
])
def test_is_palindrome(smile, dg):
    assert FreeEner.symmetry_correction(smile) == pytest.approx(dg, rel=1e-3)

def test_free_solv():
    mol_free = FreeEner('DEN', 'indene')
    mol_free.water_box_from_SMILE('C1C=Cc2ccccc12')
    mol_free.gmxsys.display()
    start_coor = pdb_manip.Coor(mol_free.gmxsys.coor_file)
    print(mol_free.gmxsys.coor_file)
    assert start_coor.num == 1625
    mol_free.equilibrate_solvent_box(em_steps=100, dt=0.002, prod_time=0.1)
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