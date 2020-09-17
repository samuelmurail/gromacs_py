#!/usr/bin/env python3

# coding: utf-8

""" Collection of function to use the antechamber toolbox.

https://ambermd.org/AmberTools.php

"""

import os
import logging
import sys

from os_command_py import os_command
from pdb_manip_py import pdb_manip

from .. import gmx

# Logging
logger = logging.getLogger(__name__)


# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

# Test folder path
MONITOR_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(MONITOR_LIB_DIR, "../test_files/"))

# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    logger.info("Ambertools cannot be found")
    REDUCE_BIN = ""
    ANTECHAMBER_BIN = ""
else:
    REDUCE_BIN = os_command.which('reduce')
    ANTECHAMBER_BIN = os_command.which('antechamber')
    ACPYPE_BIN = os_command.which('acpype')


def show_log():
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.INFO)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))


def show_debug(pdb_manip_log=True):
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.DEBUG)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))


def add_hydrogen(pdb_in, pdb_out, check_file_out=True, **reduce_options):
    """Add hydrogen to a pdb file using the ``reduce`` software:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param pdb_out: pdb output
    :type pdb_out: str

    :param check_file_out: flag to check or not if file has already been
        created. If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    :param reduce_options: Optional arguments for ``reduce``


    :Example:

    >>> from pdb_manip_py import pdb_manip
    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> add_hydrogen(pdb_in=TEST_PATH+'/phenol.pdb',\
    pdb_out=TEST_OUT+'/phenol_h.pdb') #doctest: +ELLIPSIS
    reduce -build -nuclear .../phenol.pdb
    >>> phenol_coor = pdb_manip.Coor(TEST_OUT+'/phenol_h.pdb')  \
#doctest: +ELLIPSIS
    Succeed to read file .../phenol_h.pdb ,  13 atoms found

    """
    # Check if output files exist:
    if check_file_out and os.path.isfile(pdb_out):
        logger.info("PDB files not created, {} already exist".format(
            pdb_out))
        return

    # Define reduce command:
    cmd_reduce = os_command.Command([REDUCE_BIN,
                                     "-build",
                                     "-nuclear",
                                     pdb_in], **reduce_options)

    cmd_reduce.display()
    return_code = cmd_reduce.run(display=False, out_data=True)

    filout = open(pdb_out, 'w')
    filout.write(return_code['stdout'])

    logger.info("Succeed to save file %s" % os.path.relpath(pdb_out))

    return


def antechamber(pdb_in, mol2_out, charge_model="bcc",
                check_file_out=True, **antechamber_options):
    """Compute a molecule topologie using the ``antechamber`` software:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param mol2_out: pdb output
    :type mol2_out: str

    :param charge_model: charge model
    :type charge_model: str, default="bcc"

    :param check_file_out: flag to check or not if file has already been
        created. If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    :param reduce_options: Optional arguments for ``reduce``

    Output files:

        - mol2_out
        - ANTECHAMBER_AM1BCC_PRE.AC
        - ANTECHAMBER_BOND_TYPE.AC
        - ANTECHAMBER_BOND_TYPE.AC0
        - ANTECHAMBER_AC.AC
        - ATOMTYPE.INF
        - ANTECHAMBER_AC.AC0
        - ANTECHAMBER_AM1BCC.AC

    :Example:

    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> add_hydrogen(pdb_in=TEST_PATH+'/phenol.pdb',\
    pdb_out=TEST_OUT+'/phenol_h.pdb') #doctest: +ELLIPSIS
    reduce -build -nuclear .../phenol.pdb
    >>> phenol_coor = pdb_manip.Coor(TEST_OUT+'/phenol_h.pdb')  \
#doctest: +ELLIPSIS
    Succeed to read file .../phenol_h.pdb ,  13 atoms found
    >>> antechamber(pdb_in=TEST_OUT+'/phenol_h.pdb',\
    mol2_out=TEST_OUT+'/phenol_h.mol2') #doctest: +ELLIPSIS
    antechamber -i phenol_h.pdb -fi pdb -o phenol_h.mol2 -fo \
mol2 -c bcc

    """

    # Check if output files exist:
    if check_file_out and os.path.isfile(mol2_out):
        logger.info("MOL2 files not created, {} already exist".format(
            mol2_out))
        return

    start_dir = os.path.abspath(".")
    # Create and go in out_folder:
    out_folder = os_command.get_directory(mol2_out)
    os_command.create_and_go_dir(out_folder)

    # Define reduce command:
    cmd_antechamber = os_command.Command([ANTECHAMBER_BIN,
                                          "-i", pdb_in,
                                          "-fi", str(pdb_in).split('.')[-1],
                                          "-o", mol2_out,
                                          "-fo", "mol2",
                                          "-c", charge_model],
                                         **antechamber_options)

    cmd_antechamber.display()
    cmd_antechamber.run(display=False)

    logger.info("Succeed to save file %s" % os.path.relpath(mol2_out))

    # Get absolute path:
    os.chdir(start_dir)

    return


def acpype(pdb_in, out_folder, charge_model="bcc",
           atom_type="gaff",
           check_file_out=True, **acpype_options):
    """Compute a molecule topologie using the ``antechamber`` software:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param out_folder: output folder
    :type out_folder: str

    :param charge_model: charge model
    :type charge_model: str, default="bcc"

    :param atom_type: atom_type model
    :type atom_type: str, default="gaff"

    :param check_file_out: flag to check or not if file has already been
        created. If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    :param acpype_options: Optional arguments for ``acpype``


    :Example:

    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> add_hydrogen(pdb_in=TEST_PATH+'/phenol.pdb',\
    pdb_out=TEST_OUT+'/phenol_h.pdb') #doctest: +ELLIPSIS
    reduce -build -nuclear .../phenol.pdb
    >>> phenol_coor = pdb_manip.Coor(TEST_OUT+'/phenol_h.pdb')  \
#doctest: +ELLIPSIS
    Succeed to read file .../phenol_h.pdb ,  13 atoms found

    """

    # Check if output files exist:
    if check_file_out and os.path.isfile(out_folder):
        logger.info("MOL2 files not created, {} already exist".format(
            out_folder))
        return

    # Remove tha last char if == '/'
    if out_folder[-1] == '/':
        out_folder = out_folder[:-1]

    name = out_folder.split('/')[-1]

    # Define reduce command:
    cmd_acpype = os_command.Command([ACPYPE_BIN,
                                     "-i", pdb_in,
                                     "-b", out_folder,
                                     "-c", charge_model,
                                     "-a", atom_type,
                                     "-o", "gmx"], **acpype_options)

    cmd_acpype.display()
    cmd_acpype.run(display=False)

    logger.info("Succeed to create topologie in %s" % os.path.relpath(
        out_folder))

    # Split the itp file atom_type part:
    extract_itp_atomtypes('{}.acpype/{}_GMX.itp'.format(
                          out_folder, name),
                          '{}.acpype/{}_GMX_atomtypes.itp'.format(
                          out_folder, name))

    sys_top = gmx.TopSys('{}.acpype/{}_GMX.top'.format(out_folder, name))

    fullname = '{}_GMX_atomtypes.itp'.format(name)
    include = '{}_GMX_atomtypes'.format(name)
    path = '{}.acpype/{}_GMX_atomtypes.itp'.format(out_folder, name)

    atomtypes_itp = gmx.Itp(name=include, fullname=fullname, path=path)

    sys_top.itp_list = [atomtypes_itp] + sys_top.itp_list
    sys_top.write_file('{}.acpype/{}_GMX_split.top'.format(
                                out_folder, name))

    molecule = gmx.GmxSys(name='mol',
                          coor_file='{}.acpype/{}_GMX.gro'.format(
                            out_folder, name),
                          top_file='{}.acpype/{}_GMX_split.top'.format(
                            out_folder, name))

    return molecule


def make_amber_top_mol(pdb_in, res_name, charge_model="bcc",
                       atom_type="gaff"):

    full_coor = pdb_manip.Coor(pdb_in)
    # Select coor:
    mol_coor = full_coor.select_part_dict(selec_dict={'res_name': [res_name]})
    mol_coor.write_pdb(res_name+'.pdb')

    # Add hydrogens:
    add_hydrogen(res_name+'.pdb', res_name+'_h.pdb')

    # Get only one molecule
    mol_h_coor = pdb_manip.Coor(res_name+'_h.pdb')
    res_list = mol_h_coor.get_attribute_selection(attribute='uniq_resid')
    mol_uniq_coor = mol_h_coor.select_part_dict(
        selec_dict={'uniq_resid': [res_list[0]]})
    # Save coordinates:
    mol_uniq_coor.write_pdb(res_name+'_h_unique.pdb')

    # Compute topologie with acpype:
    gmxsys = acpype(res_name+'_h_unique.pdb', res_name, charge_model="bcc",
                    atom_type="gaff")

    return({'GmxSys': gmxsys,
            'coor': os_command.full_path_and_check(res_name+'_h.pdb'),
            'num': len(res_list)})


def extract_itp_atomtypes(itp_in, itp_atomtypes_out):

    field = None
    atom_types_list = "[ atomtypes ]\n"

    with open(itp_in) as itpfile:
        for line in itpfile:

            if line.startswith("["):
                # Remove space and [ ]
                field = line.strip()[1:-1].strip()
                # print(field)
                continue
            if field == 'atomtypes':
                atom_types_list += line

    filout = open(itp_atomtypes_out, "w")
    filout.write(atom_types_list)
    filout.close()

    # Read and write itp to remove the [ atomtypes ] part
    fullname = (itp_in.split("/")[-1])
    include = fullname.split(".")[0]
    path = os_command.full_path_and_check(itp_in)

    top_mol = gmx.Itp(name=include, fullname=fullname, path=path)
    top_mol.write_file(path)
