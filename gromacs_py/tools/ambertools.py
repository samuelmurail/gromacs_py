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
    # As the mabertools module is not a hard dependency
    # gromacs_py should be loaded even if ambertools 
    # is note installes
    try:
        REDUCE_BIN = os_command.which('reduce')
        ANTECHAMBER_BIN = os_command.which('antechamber')
        ACPYPE_BIN = os_command.which('acpype', 'acpype.py')
    except IOError:
        logger.warning('Reduce or antechamber or acpype binary cannot be'
                       ' found. Install it using conda \n conda install '
                       'acpype')


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


def smile_to_pdb(smile, pdb_out, mol_name,
                 method_3d='rdkit', iter_num=5000):
    """
    """

    if method_3d == 'openbabel':

        from openbabel import pybel

        conf = pybel.readstring("smi", smile)
        # Get charge
        charge = conf.charge
        conf.make3D(forcefield='mmff94', steps=iter_num)
        conf.write(format='pdb', filename=pdb_out, overwrite=True)

    elif method_3d == 'rdkit':

        from rdkit.Chem import AllChem as Chem

        conf = Chem.MolFromSmiles(smile)
        conf = Chem.AddHs(conf)
        # Get charge
        charge = Chem.GetFormalCharge(conf)
        Chem.EmbedMolecule(conf)
        Chem.MMFFOptimizeMolecule(
            conf, mmffVariant='MMFF94', maxIters=iter_num)
        Chem.MolToPDBFile(conf, filename=pdb_out)

    # Change resname of pdb file to `self.mol_name`
    coor = pdb_manip.Coor(pdb_out)
    index_list = coor.get_index_selection(selec_dict={'res_name': ['UNL']})
    coor.change_index_pdb_field(index_list, change_dict={
        'res_name': mol_name})
    coor.write_pdb(pdb_out, check_file_out=False)

    return(charge)


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


def add_hydrogen_rdkit(pdb_in, smile, pdb_out):
    """Add hydrogen to a pdb file using the ``rdkit`` library:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param pdb_out: pdb output
    :type pdb_out: str

    :Example:

    >>> from pdb_manip_py import pdb_manip
    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> # print(TEST_OUT)
    >>> add_hydrogen_rdkit(pdb_in=os.path.join(TEST_PATH,'four_phenol.pdb'),\
smile="C1=CC=C(C=C1)O",\
pdb_out=os.path.join(TEST_OUT,'four_phenol_h.pdb')) #doctest: +ELLIPSIS
    Succeed to read file ...four_phenol.pdb ,  28 atoms found
    Succeed to save file ...four_phenol_0.pdb
    Succeed to read file ...four_phenol_0.pdb ,  7 atoms found
    Succeed to read file ...four_phenol_0_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_0_h.pdb
    Succeed to read file ...four_phenol_0_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_0_h.pdb
    Succeed to save file ...four_phenol_1.pdb
    Succeed to read file ...four_phenol_1.pdb ,  7 atoms found
    Succeed to read file ...four_phenol_1_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_1_h.pdb
    Succeed to read file ...four_phenol_1_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_1_h.pdb
    Succeed to save file ...four_phenol_2.pdb
    Succeed to read file ...four_phenol_2.pdb ,  7 atoms found
    Succeed to read file ...four_phenol_2_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_2_h.pdb
    Succeed to read file ...four_phenol_2_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_2_h.pdb
    Succeed to save file ...four_phenol_3.pdb
    Succeed to read file ...four_phenol_3.pdb ,  7 atoms found
    Succeed to read file ...four_phenol_3_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_3_h.pdb
    Succeed to read file ...four_phenol_3_h.pdb ,  13 atoms found
    Succeed to save file ...four_phenol_3_h.pdb
    Succeed to save concat file:  ...four_phenol_h.pdb
    0
    >>> phenol_coor = pdb_manip.Coor(os.path.join(TEST_OUT,'four_phenol_h.pdb'))\
#doctest: +ELLIPSIS
    Succeed to read file ...four_phenol_h.pdb ,  52 atoms found

    """

    full_coor = pdb_manip.Coor(pdb_in)

    # Change first residue to 1, as it is the res from rdkit output
    res_list = full_coor.get_attribute_selection(attribute='res_num')
    pdb_list = []
    pdb_del_list = []

    for i, res in enumerate(res_list):
        # Extract residue
        one_res = full_coor.select_part_dict(
            selec_dict={'res_num': [res]})
        out_pdb = '{}_{}.pdb'.format(pdb_in[:-4], i)
        one_res.write_pdb(out_pdb)
        pdb_del_list.append(out_pdb)

        # Add hydrogens
        out_h_pdb = '{}_{}_h.pdb'.format(pdb_in[:-4], i)
        charge = add_hydrogen_rdkit_one_mol(out_pdb, smile, out_h_pdb)
        pdb_list.append(out_h_pdb)
        pdb_del_list.append(out_h_pdb)

        # Update residue number:
        one_res_h = pdb_manip.Coor(out_h_pdb)
        one_res_h.change_pdb_field(
            change_dict={"res_num": i + 1})
        one_res_h.write_pdb(out_h_pdb, check_file_out=False)


    pdb_manip.Coor.concat_pdb(*pdb_list,
                              pdb_out = pdb_out)

    # Delete intermediate files
    for pdb in pdb_del_list:
        os_command.delete_file(pdb)

    return charge

def add_hydrogen_rdkit_one_mol(pdb_in, smile, pdb_out):
    """Add hydrogen to a pdb file using the ``rdkit`` library:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param pdb_out: pdb output
    :type pdb_out: str

    :Example:

    >>> from pdb_manip_py import pdb_manip
    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> # print(TEST_OUT)
    >>> add_hydrogen_rdkit_one_mol(pdb_in=os.path.join(TEST_PATH,'phenol.pdb'),\
smile="C1=CC=C(C=C1)O",\
pdb_out=os.path.join(TEST_OUT,'phenol_h.pdb')) #doctest: +ELLIPSIS
    Succeed to read file ...phenol_h.pdb ,  13 atoms found
    Succeed to save file ...phenol_h.pdb
    0
    >>> phenol_coor = pdb_manip.Coor(os.path.join(TEST_OUT,'phenol_h.pdb'))\
#doctest: +ELLIPSIS
    Succeed to read file .../phenol_h.pdb ,  13 atoms found

    .. warning:

        

    """

    try:
        from rdkit.Chem import AllChem as Chem
        # from rdkit.Chem import rdMolAlign as align
    except ImportError:
        logger.error('Could not load rdkit \nInstall it using conda:\n'
                     'conda install -c conda-forge rdkit')
        sys.exit(1)

    lig_pdb = Chem.MolFromPDBFile(pdb_in, removeHs=True)
    lig_smile = Chem.MolFromSmiles(smile)

    # Need to count the number of molecule
    pdb_atom_num = lig_pdb.GetNumAtoms()
    smile_atom_num = lig_smile.GetNumAtoms()

    # If more than one molecule, add them
    # in the smile string
    if pdb_atom_num != smile_atom_num:
        mol_num = pdb_atom_num/smile_atom_num
        smile_list = [smile] * int(mol_num)
        lig_smile = Chem.MolFromSmiles('.'.join(smile_list))

    # Assign bond order on pdb using smile informations
    newMol = Chem.AssignBondOrdersFromTemplate(lig_smile, lig_pdb)

    # Add hydrogens to lig_smile, using lig_pdb as coordinates constraints
    # This way is better than adding h to lig_pdb
    # because acpype experience some issue
    lig_smile_h = Chem.AddHs(lig_smile)
    Chem.ConstrainedEmbed(lig_smile_h, newMol)
    Chem.MolToPDBFile(lig_smile_h, pdb_out)

    # OLD WAY :
    # Add hydrogens
    # newMol_h = Chem.AddHs(newMol)

    # # Need to define how to match atoms form pdb to smile
    # match_atom = newMol_h.GetSubstructMatch(lig_smile)
    # cmap = {match_atom[i]: lig_pdb.GetConformer().
    #         GetAtomPosition(match_atom[i]) for i in range(len(match_atom))}

    # # Hydrogens coordinates need to be computed
    # # While keeping heavy atoms coordinates
    # Chem.EmbedMolecule(newMol_h, coordMap=cmap)

    # # Align new coordinates to old one
    # align.AlignMol(newMol_h, newMol,
    #                atomMap=[[i, i] for i in range(len(match_atom))])
    #  Save coordinates
    # Chem.MolToPDBFile(newMol_h, pdb_out)

    # Change UNL residue name to original one
    coor_start = pdb_manip.Coor(pdb_in)
    res_name_list_start = coor_start.get_attribute_selection(
        attribute='res_name')
    if len(res_name_list_start) > 1:
        res_name_list_start.remove('UNL')

    coor = pdb_manip.Coor(pdb_out)
    index_list = coor.get_index_selection(selec_dict={'res_name': ['UNL']})
    coor.change_index_pdb_field(index_list, change_dict={
        'res_name': res_name_list_start[0]})
    coor.write_pdb(pdb_out, check_file_out=False)

    # Return charge
    return Chem.GetFormalCharge(lig_smile)

    # Need to fix the residue number


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
           atom_type="gaff", net_charge=None,
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

    cmd_list = [ACPYPE_BIN,
                "-i", pdb_in,
                "-b", out_folder,
                "-c", charge_model,
                "-a", atom_type,
                "-o", "gmx"]

    if net_charge is not None:
        cmd_list += ["-n", str(net_charge)]

    # Define reduce command:
    cmd_acpype = os_command.Command(cmd_list, **acpype_options)

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
    # Add position restraints
    # Get heavy atoms name:
    atom_list = []
    for key, value in sys_top.itp_list[-1].top_mol_list[0].atom_dict.items():
        if not value['atom_name'].startswith('H'):
            atom_list.append(value['atom_name'])

    sys_top.add_posre(posre_name="POSRES", selec_dict={
            'atom_name': atom_list},
            fc=[1000, 1000, 1000])
    sys_top.add_posre(posre_name="POSRES_LIG", selec_dict={
            'atom_name': atom_list},
            fc=[1000, 1000, 1000])

    sys_top.add_posre(posre_name="POSRES_CA", selec_dict={
            'atom_name': atom_list},
            fc=[1000, 1000, 1000])

    sys_top.add_posre(posre_name="POSRES_CA_LOW", selec_dict={
            'atom_name': atom_list},
            fc=[100, 100, 100])

    sys_top.write_file('{}.acpype/{}_GMX_split.top'.format(
                                out_folder, name))

    molecule = gmx.GmxSys(name='mol',
                          coor_file='{}.acpype/{}_GMX.gro'.format(
                            out_folder, name),
                          top_file='{}.acpype/{}_GMX_split.top'.format(
                            out_folder, name))

    return molecule


def make_amber_top_mol(pdb_in, res_name, charge, charge_model="bcc",
                       atom_type="gaff", remove_h=True):

    full_coor = pdb_manip.Coor(pdb_in)
    # Select coor:
    mol_coor = full_coor.select_part_dict(selec_dict={'res_name': [res_name]})

    # Remove hydrogens:
    if remove_h:
        to_del_list = []
        for atom_num, atom in mol_coor.atom_dict.items():
            if atom['name'][0] == 'H':
                to_del_list.append(atom_num)
        mol_coor.del_atom_index(to_del_list)

    mol_coor.write_pdb(res_name + '.pdb')

    # Add hydrogens:
    add_hydrogen(res_name + '.pdb', res_name + '_h.pdb')

    # Get only one molecule
    mol_h_coor = pdb_manip.Coor(res_name + '_h.pdb')
    res_list = mol_h_coor.get_attribute_selection(attribute='uniq_resid')
    mol_uniq_coor = mol_h_coor.select_part_dict(
        selec_dict={'uniq_resid': [res_list[0]]})
    # Save coordinates:
    mol_uniq_coor.write_pdb(res_name + '_h_unique.pdb')

    # Compute topologie with acpype:
    gmxsys = acpype(res_name + '_h_unique.pdb', res_name,
                    net_charge=charge,
                    charge_model="bcc",
                    atom_type="gaff")

    return({'GmxSys': gmxsys,
            'coor': os_command.full_path_and_check(res_name + '_h.pdb'),
            'num': len(res_list)})


def make_amber_top_mol_rdkit(pdb_in, res_name, smile, charge_model="bcc",
                             atom_type="gaff", remove_h=True):

    full_coor = pdb_manip.Coor(pdb_in)
    # Select coor:
    mol_coor = full_coor.select_part_dict(selec_dict={'res_name': [res_name]})

    # Change first residue to 1, as it is the res from rdkit output
    res_list = mol_coor.get_attribute_selection(attribute='res_num')
    index_all_list = []
    for res in res_list:
        index_list = mol_coor.get_index_selection(
            selec_dict={'res_num': [res]})
        index_all_list.append(index_list)

    for i, index_list in enumerate(index_all_list):
        mol_coor.change_index_pdb_field(index_list, change_dict={'res_num': i+1})

    #index_list = mol_coor.get_index_selection(
    #    selec_dict={'res_num': [res_list[0]]})
    #mol_coor.change_index_pdb_field(index_list, change_dict={'res_num': 1})

    # Remove hydrogens:
    if remove_h:
        to_del_list = []
        for atom_num, atom in mol_coor.atom_dict.items():
            if atom['name'].startswith('H'):
                to_del_list.append(atom_num)
        mol_coor.del_atom_index(to_del_list)

    mol_coor.write_pdb(res_name + '.pdb', check_file_out=False)

    # Add hydrogens:
    # add_hydrogen(res_name+'.pdb', res_name+'_h.pdb')
    charge = add_hydrogen_rdkit(res_name + '.pdb',
                                smile,
                                res_name + '_h.pdb')

    # Get only one molecule
    mol_h_coor = pdb_manip.Coor(res_name + '_h.pdb')
    res_list = mol_h_coor.get_attribute_selection(attribute='uniq_resid')
    mol_uniq_coor = mol_h_coor.select_part_dict(
        selec_dict={'uniq_resid': [res_list[0]]})
    # Save coordinates:
    mol_uniq_coor.write_pdb(res_name + '_h_unique.pdb', check_file_out=False)

    # Compute topologie with acpype:
    gmxsys = acpype(res_name+'_h_unique.pdb', res_name,
                    net_charge=charge,
                    charge_model="bcc",
                    atom_type="gaff")

    return({'GmxSys': gmxsys,
            'coor': os_command.full_path_and_check(res_name + '_h.pdb'),
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
