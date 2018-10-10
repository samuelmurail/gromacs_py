#!/usr/bin/env python3
# coding: utf-8
###################################
#########     PDB2PQR    ##########
###################################

__author__ = "Samuel Murail"

import sys
import  os
# Needed for doctest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import tools.osCommand as osCommand
import tools.pdb_manip as pdb_manip

#PDB2PQR_BIN="/Users/murail/Documents/test/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py"
#PDB2PQR_BIN="/Users/murail/Documents/Software/pdb2pqr/pdb2pqr.py"
#PDB2PQR_BIN=os.path.expanduser("~/Documents/git/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py")

PDB2PQR_MOD_DIRNAME = os.path.dirname(os.path.abspath(__file__))
#FORCEFIELD_PATH=os.path.abspath(GROMACS_MOD_DIRNAME+"/template/")
#PDB2PQR_BIN=os.path.expanduser(PDB2PQR_MOD_DIRNAME+"/../../apbs-pdb2pqr/pdb2pqr/pdb2pqr.py")

PDB2PQR_BIN='pdb2pqr.py'


def compute_pdb2pqr(pdb_in, pdb_out, ff = "CHARMM", check_file_out = True):
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
    >>> pdb2pqr.compute_pdb2pqr('../test/input/4n1m.pdb', '../test/output/pdb2pqr_test/4n1m.pqr')
    Succeed to read file ../test/input/4n1m.pdb ,  2530 atoms found
    Succeed to save file ../test/output/pdb2pqr_test/tmp_pdb2pqr.pdb
    pdb2pqr.py --ff CHARMM --ffout CHARMM --chain ../test/output/pdb2pqr_test/tmp_pdb2pqr.pdb ../test/output/pdb2pqr_test/4n1m.pqr
    0
    >>> prot_coor = pdb_manip.coor()
    >>> prot_coor.read_pdb('../test/output/pdb2pqr_test/4n1m.pqr', pqr_format = True)
    Succeed to read file ../test/output/pdb2pqr_test/4n1m.pqr ,  2548 atoms found
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
    if check_file_out and osCommand.check_file_and_create_path(pdb_out):
        #print("pdb2pqr not launched",pdb_out,"already exist")
        return(pdb_out)

    out_folder = osCommand.get_directory(pdb_out)
    #print("out_folder", out_folder)

    # WARING :
    # Many bugs are due to the REMARK field in pdb2pqr
    # The 2 following steps remove the REMARK field of the pdb

    tmp_coor = pdb_manip.coor()
    tmp_coor.read_pdb(pdb_in)

    # Remove HETATM
    no_hetatm_pdb = tmp_coor.select_part_dict({'field':'ATOM'})
    no_hetatm_pdb.write_pdb(out_folder+"/tmp_pdb2pqr.pdb")


    cmd_pdb2pqr = osCommand.command([PDB2PQR_BIN, 
        "--ff", ff, 
        "--ffout", ff,
        "--chain",
        out_folder+"/tmp_pdb2pqr.pdb", pdb_out])

    cmd_pdb2pqr.display()
    out_data = cmd_pdb2pqr.run()
    osCommand.delete_file(out_folder+"/tmp_pdb2pqr.pdb")

    return(out_data)

if __name__ == "__main__":

    import doctest
    import shutil
    doctest.testmod()
    # Erase all test files
    shutil.rmtree('../test/output/pdb2pqr_test', ignore_errors=True)


