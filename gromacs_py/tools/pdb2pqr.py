#!/usr/bin/env python3
# coding: utf-8
"""
###################################
#########     PDB2PQR    ##########
###################################
"""

__author__ = "Samuel Murail"

import  os
# Needed for doctest
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import tools.os_command as os_command
import tools.pdb_manip as pdb_manip

PDB2PQR_MOD_DIRNAME = os.path.dirname(os.path.abspath(__file__))
PDB2PQR_BIN = 'pdb2pqr.py'

# Test folder path
PQR_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(PQR_LIB_DIR+"/../test/input/")
TEST_OUT = 'gromacs_py_test_out/pdb2pqr_test'


def compute_pdb2pqr(pdb_in, pdb_out, ff="CHARMM", check_file_out=True):
    """
    Use pdb2pqr to define protonation state of each residue of a protein.

    :param pdb_in: path of input pdb file
    :type pdb_in: str

    :param pdb_out: path of output pdb file
    :type pdb_out: str

    :param ff: forcefield nomenclature for atom names
    :type ff: str, optional, default="CHARMM"

    :param check_file_out: flag to check or not if file has already been created.
        If the file is present then the command break.
    :type check_file_out: bool, optional, default=True


    :Example:

    >>> import tools.pdb_manip as pdb_manip
    >>> import tools.pdb2pqr as pdb2pqr
    >>> # Compute protonation with pdb2pqr:
    >>> pdb2pqr.compute_pdb2pqr(TEST_PATH+'/4n1m.pdb', TEST_OUT+'/4n1m.pqr') #doctest: +ELLIPSIS
    Succeed to read file ...test/input/4n1m.pdb ,  2530 atoms found
    Succeed to save file gromacs_py_test_out/pdb2pqr_test/tmp_pdb2pqr.pdb
    pdb2pqr.py --ff CHARMM --ffout CHARMM --chain gromacs_py_test_out/pdb2pqr_test/tmp_pdb2pqr.pdb gromacs_py_test_out/pdb2pqr_test/4n1m.pqr
    0
    >>> prot_coor = pdb_manip.Coor()
    >>> prot_coor.read_pdb(TEST_OUT+'/4n1m.pqr', pqr_format = True)
    Succeed to read file gromacs_py_test_out/pdb2pqr_test/4n1m.pqr ,  2548 atoms found
    >>> HSD_index = prot_coor.get_index_selection({'res_name' : ['HSD'], 'name':['CA']})
    >>> print(len(HSD_index))
    5
    >>> HSE_index = prot_coor.get_index_selection({'res_name' : ['HSE'], 'name':['CA']})
    >>> print(len(HSE_index))
    0
    >>> HSP_index = prot_coor.get_index_selection({'res_name' : ['HSP'], 'name':['CA']})
    >>> print(len(HSP_index))
    0

    .. note::
        Idealy I would need a pdb file with 3 different histidine protonation. I couldn't find one.

    """

    #print("Compute pdb2pqr on",pdb_in)

    # Check if output files exist and create directory:
    if check_file_out and os_command.check_file_and_create_path(pdb_out):
        #print("pdb2pqr not launched",pdb_out,"already exist")
        return pdb_out

    out_folder = os_command.get_directory(pdb_out)
    #print("out_folder", out_folder)

    # WARING :
    # Many bugs are due to the REMARK field in pdb2pqr
    # The 2 following steps remove the REMARK field of the pdb

    tmp_coor = pdb_manip.Coor()
    tmp_coor.read_pdb(pdb_in)

    # Remove HETATM
    no_hetatm_pdb = tmp_coor.select_part_dict({'field':'ATOM'})
    no_hetatm_pdb.write_pdb(out_folder+"/tmp_pdb2pqr.pdb")


    cmd_pdb2pqr = os_command.Command([PDB2PQR_BIN,
                                      "--ff", ff,
                                      "--ffout", ff,
                                      "--chain",
                                      out_folder+"/tmp_pdb2pqr.pdb", pdb_out])

    cmd_pdb2pqr.display()
    out_data = cmd_pdb2pqr.run()
    os_command.delete_file(out_folder+"/tmp_pdb2pqr.pdb")

    return out_data

if __name__ == "__main__":

    import doctest
    import shutil
    doctest.testmod()
    # Erase all test files
    shutil.rmtree('../test/output/pdb2pqr_test', ignore_errors=True)
