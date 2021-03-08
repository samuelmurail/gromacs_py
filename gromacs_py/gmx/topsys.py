#!/usr/bin/env python3
# coding: utf-8
##################################
# #######   GROMACS 5   ##########
##################################


import os
import logging
import sys

from shutil import copy as shutil_copy

from os_command_py import os_command

from .itp import Itp

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

# Logging
logger = logging.getLogger(__name__)


# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    logger.info("Gromacs cannot be found")
    GMX_BIN = ""
    gmx_version = ""
else:
    GMX_BIN = os_command.which('gmx')
    gmx_version = os_command.get_gmx_version()
    logger.info("Gromacs version is {}".format(gmx_version))

# Get local gmx path
GMX_PATH = "/".join(GMX_BIN.split("/")[:-2])

# Get library path
GROMACS_MOD_DIRNAME = os.path.dirname(os.path.abspath(__file__))

FORCEFIELD_PATH_LIST = [os.path.join(GROMACS_MOD_DIRNAME, "template")]
# GMXLIB_var should not include GMXPATH top dir, it could induce some conflict
GMXLIB_var = os.path.join(GROMACS_MOD_DIRNAME, "template")

# Check if GMXLIB env variable is defines if yes add it to forcefield path
if 'GMXLIB' in os.environ:
    # Needed for pytest in monitor.py, otherwise add twice GROMACS_MOD_DIRNAME
    if os.environ['GMXLIB'] not in FORCEFIELD_PATH_LIST:
        FORCEFIELD_PATH_LIST.append(os.environ['GMXLIB'])
        GMXLIB_var += ":" + os.environ['GMXLIB']

FORCEFIELD_PATH_LIST.append(os.path.join(GMX_PATH, "share/gromacs/top"))


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

    # Show log of top sys:
    from . import itp
    itp.show_log()


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


class TopSys:
    """Topologie base on gromacs .top :
    #include forcefield

    All name and full path of itp are save
    [ system ] -> Name
    [ molecules ] -> Composition

    :param path: topologie file path
    :type path: str

    :param forcefield: name of the focefield
    :type forcefield: str

    :param itp_list: list of the itp object
    :type itp_list: list

    :param mol_comp: molecular composition
    :type mol_comp: list

    :param name: name of the system
    :type tpr: str

    :param folder: path of the top file folder
    :type folder: str

    :param include_itp: Flag indicating if the topologie include a molecule
        topologie
    :type include_itp: bool
    """

    def __init__(self, top_in):

        self.path = os_command.full_path_and_check(top_in)
        self.forcefield = None
        self.itp_list = []
        self.mol_comp = []
        self.name = ""
        self.folder = os_command.get_directory(top_in)
        self.include_itp = False
        # Intermolecular restraints (only for Free ener):
        self.inter_interact = False
        self.inter_bond_list = []
        self.inter_angl_list = []
        self.inter_dihe_list = []

        self.read_file(top_in)
        if self.include_itp:
            logger.info("Rewrite topologie: {}".format(top_in))
            self.write_file(top_in)

    def read_file(self, top_in):

        field = None
        ifdef = False
        bond_list = []
        angle_list = []
        dihed_list = []

        with open(top_in) as topfile:
            for line in topfile:
                # print("line: ",line)
                # Check #ifdef field:
                if line[:6] == '#ifdef':
                    ifdef = True
                if line[:6] == '#endif':
                    ifdef = False
                # Check include files
                if not ifdef and line.startswith('#include'):
                    # get the file name to include:
                    file_name = line.split()[1][1:-1]
                    # remove '"' , path and .itp:
                    include = (file_name.split("/")[-1]).split(".")[0]

                    # Look if the itp is in the top path:
                    # print("Path 1: ",self.folder+"/"+file_name)
                    # print("Path 2: ",FORCEFIELD_PATH_LIST+"/"+file_name)

                    itp_found = False
                    if os_command.check_file_exist(
                            os.path.join(self.folder, file_name)):
                        path = os.path.abspath(os.path.join(self.folder,
                                                            file_name))
                        itp_found = True
                    else:
                        for forcefield in FORCEFIELD_PATH_LIST:
                            if os_command.check_file_exist(
                                    os.path.join(forcefield, file_name)):
                                path = os.path.abspath(
                                    os.path.join(forcefield, file_name))
                                itp_found = True
                                break
                    if not itp_found:
                        raise IOError('Itp ' + file_name + ' not found')

                    # print("name =", include, "fullname =", file_name,
                    # "path =",path)
                    if include == "forcefield":
                        self.forcefield = {'name': file_name.split('.')[0],
                                           'fullname': file_name,
                                           'path': path}
                    else:
                        self.add_mol_itp(path)
                # Check field
                elif not ifdef and line.startswith("["):
                    # Remove space and [ ]
                    field = line.strip()[1:-1].strip()
                    # print(field)
                    # As intermolecular_inter field is empty
                    # (filled with bonds, angles, ...)
                    # It has to be checked now:
                    if field == 'intermolecular_interactions':
                        self.inter_interact = True
                    continue
                # Check in the field :
                elif (not ifdef and not line[0].startswith(";") and
                      line.strip() != ""):
                    if field == 'moleculetype':
                        # in the case where mol param are present in the top
                        # file, convert the top to an itp
                        # If mol_name variable is defined, give to the mol
                        # param this name
                        name_itp = os.path.basename(top_in)[:-3] + "itp"

                        logger.info("Molecule topologie present in {} "
                                    ", extract the topologie in a separate"
                                    " file: {}".format(top_in, name_itp))
                        # Store and write the itp file:
                        top_itp = Itp(name=name_itp, fullname=name_itp,
                                      path=os.path.abspath(top_in))
                        top_itp.write_file(os.path.join(self.folder, name_itp))
                        top_itp.display()
                        # Add the itp to the itp_list
                        self.add_mol_itp(os.path.join(self.folder, name_itp))
                        self.include_itp = True
                    # Name of the system
                    elif field == 'system':
                        self.name = line.strip()
                    # Molecule composition of the system
                    elif field == 'molecules':
                        line_list = line.strip().split()
                        self.mol_comp.append({'name': line_list[0],
                                              'num': line_list[1]})
                    # Following command concern
                    # intermolecular interactions
                    elif field == 'bonds':
                        bond_list.append(line.strip().split())
                    elif field == 'angles':
                        angle_list.append(line.strip().split())
                    elif field == 'dihedrals':
                        dihed_list.append(line.strip().split())
        if self.inter_interact:
            self.add_intermolecular_restr(bond_list=bond_list,
                                          angle_list=angle_list,
                                          dihed_list=dihed_list)

    def display(self):
        if self.forcefield:
            print("Forcefield include :\n", self.forcefield['name'])
        else:
            print('No forcefield defined !')
        for itp in self.itp_list:
            itp.display()
        print("Mol List:")
        for mol in self.mol_comp:
            print("   * {} {}".format(mol['num'], mol['name']))

        print("Mol Name:\n", self.name)

    def write_file(self, top_out):
        filout = open(top_out, 'w')

        filout.write("; Topologie file created by " + __name__ + "\n")
        if self.forcefield:
            filout.write("; Forcefield: \n")
            filout.write("#include \"" + self.forcefield['fullname'] + "\"\n")
        # print include files
        filout.write("\n; Itp to include: \n")
        for itp in self.itp_list:
            filout.write("#include \"" + itp.fullname + "\"\n")
            # Check if the include is in the ff folder:
            # if self.forcefield['fullname'].split("/")[0] == \
            #       itp.fullname.split("/")[0]:
            #    filout.write("#include \""+itp.fullname+"\"\n")
            # else:
            #    filout.write("#include \""+itp.fullname+"\"\n")
        filout.write("\n[ system ]\n; Name\n")
        filout.write(self.name + "\n")
        filout.write("\n[ molecules ]\n; Compound        #mols\n")
        for mol in self.mol_comp:
            filout.write(mol['name'] + " \t" + mol['num'] + "\n")
        filout.write("\n")
        # INTER
        if self.inter_interact:
            filout.write("[ intermolecular_interactions ]\n")
            if self.inter_bond_list:
                filout.write("\n[ bonds ]\n;  ai    aj    type   "
                             " bA    kA    bB    kB\n")
                for param in self.inter_bond_list:
                    filout.write(
                        "{:>6}{:>6}{:>6}{:>13}{:>13}{:>13}{:>13}\n".format(
                            param['ai'], param['aj'], param['type'],
                            param['rA'], param['kA'],
                            param['rB'], param['kB']))
            if self.inter_angl_list:
                filout.write("\n[ angles ]\n;  ai    aj    ak  "
                             "  type    thA    kA    thB    kB\n")
                for param in self.inter_angl_list:
                    filout.write("{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}{:>13}"
                                 "{:>13}\n".format(
                                    param['ai'], param['aj'],
                                    param['ak'], param['type'],
                                    param['thA'], param['kA'],
                                    param['thB'], param['kB']))
            if self.inter_dihe_list:
                filout.write("\n[ dihedrals ]\n;  ai    aj    ak    al"
                             "    type    phiA    kA    phiB    kB\n")
                for param in self.inter_dihe_list:
                    filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}"
                                 "{:>13}{:>13}\n".format(
                                    param['ai'], param['aj'],
                                    param['ak'], param['al'],
                                    param['type'],
                                    param['thA'], param['kA'],
                                    param['thB'], param['kB']))

        filout.close()
        self.copy_dependancies(os_command.get_directory(top_out))

    def charge(self):
        """Get the charge of the system
        """

        self._charge = 0
        # Look for mol name
        for mol in self.mol_comp:
            name = mol['name']
            num = int(mol['num'])
            # print("Mol top info: ",name,num)
            # Search mol name in itp:
            for local_itp in self.itp_list:
                itp_charge = local_itp.charge(name)
                if itp_charge is not None:
                    logger.debug("Get charge of {} : {} total charge:"
                                 " {}".format(name, itp_charge,
                                              itp_charge * num))
                    self._charge += num * itp_charge
                    break
        return self._charge

    def prot_res_num(self, selection="Protein"):
        """Compute the residue number of a selection
        """

        self.res_num = 0
        # Look for mol name
        for mol in self.mol_comp:
            name = mol['name']
            num = int(mol['num'])
            logger.info("{} : {}".format(name, num))
            if name.find(selection) != -1:
                # Search mol name in itp:
                for itp in self.itp_list:
                    # print(itp.name)
                    itp_res_num = itp.res_num(name)
                    if itp_res_num is not None:
                        logger.info("Get Res num of {} : {}\n"
                                    "Total number of residue: {}".format(
                                        name, itp_res_num, num * itp_res_num))
                        self.res_num += num * itp_res_num
                        break
        return self.res_num

    def mol_num(self, name):
        """Get the number of the molecule "name"
        """

        mol_num = 0

        for mol in self.mol_comp:
            if mol['name'] == name:
                mol_num = mol_num + int(mol['num'])
        return mol_num

    def add_posre(self, posre_name="POSRE_CA",
                  selec_dict={'atom_name': ['CA']},
                  fc=[1000, 1000, 1000]):
        """Add position restraint based on the selection for each itp
        """
        for mol in self.mol_comp:
            for itp in self.itp_list:
                itp.add_posre(mol_name=mol['name'], posre_name=posre_name,
                              selec_dict=selec_dict, fc=fc)

    def add_mol_itp(self, mol_itp_file):
        """Add a molecule itp in the topologie itp_list.
        """
        fullname = (mol_itp_file.split("/")[-1])
        include = fullname.split(".")[0]
        path = os_command.full_path_and_check(mol_itp_file)
        mol_itp = Itp(name=include, fullname=fullname, path=path)

        # 1. Get new mol name:
        mol_name_list = []
        for top_mol in mol_itp.top_mol_list:
            mol_name_list.append(top_mol.name)

        # 2. Check if it already present
        present = False
        for itp in self.itp_list:
            for top_mol in itp.top_mol_list:
                if top_mol.name in mol_name_list:
                    present = True
                    break
            if present:
                break

        # Add the itp if not present
        if not present:
            self.itp_list.append(mol_itp)

    def add_mol(self, mol_name, mol_itp_file, mol_num):
        """Add a molecule in the topologie (composition and itp_list)
        """

        logger.info("Add {} mol {}".format(mol_num, mol_itp_file))
        self.mol_comp.append({'name': mol_name, 'num': str(mol_num)})
        self.add_mol_itp(mol_itp_file)

    def remove_ion(self, ion_name_list):
        """Remove a molecule from the topologie (composition and itp_list)
        """

        logger.info("Remove mol(s) : {}".format(' '.join(ion_name_list)))
        itp_remove_list = []
        mol_remove_list = []

        for itp in self.itp_list:
            for mol in itp.top_mol_list:
                for atom in mol.atom_dict.values():
                    if atom['atom_name'] in ion_name_list and (
                            len(itp.top_mol_list) == 1):
                        # print(itp.name, 'YOYOYO')
                        itp_remove_list.append(itp)
                        mol_remove_list.append(mol.name)

        for itp in itp_remove_list:
            self.itp_list.remove(itp)

        mol_comp = []
        for mol in self.mol_comp:
            if mol['name'] not in mol_remove_list:
                mol_comp.append(mol)

        self.mol_comp = mol_comp

    def add_atomtypes(self, new_atomtypes):
        """ Add atomtypes in a topologie.

        :param new_atomtypes: path of the atomtype itp file
        :type new_atomtypes: str

        """

        # check if and atomtypes file exists:
        # self.display()
        atom_type = False
        for itp_file in self.get_include_no_posre_file_list():
            if itp_file.split("/")[-1].endswith('atomtypes.itp'):
                # print('atomtypes file present', itp_file)
                atom_type = True
                atomtype_path = itp_file
                break

        # If it doesn't exist create it
        if not atom_type:
            fullname = (new_atomtypes.split("/")[-1])
            include = fullname.split(".")[0]
            path = os_command.full_path_and_check(new_atomtypes)
            atomtype_itp = Itp(name=include, fullname=fullname, path=path)

            self.itp_list = [atomtype_itp] + self.itp_list
        # If it does exist add the new atomtypes in the first one:
        else:
            # First extract old atom types:
            field = None
            atom_dict = {}
            name_list = []
            with open(atomtype_path) as file:
                for line in file:
                    if line.strip().startswith("["):
                        # Remove space and [ ], remove also comments
                        field = line.replace(" ", "").split("]")[0][1:]
                        continue

                    if (line[0] != ";" and line[0] != "#" and
                            line.strip() != ""):
                        # Remove commentary in the line
                        line_comment = line.split(';')
                        line_list = line_comment[0].split()

                        if field == 'atomtypes':
                            name_list.append(line_list[0])
                            name, bond_type, mass, charge, ptype, sigma, \
                                epsilon = line_list[:7]
                            atom_dict[name] = {'bond_type': bond_type,
                                               'mass': float(mass),
                                               'charge': float(charge),
                                               'ptype': ptype,
                                               'sigma': float(sigma),
                                               'epsilon': float(epsilon)}
                            name_list.append(name)

            # Second extract new atom types and
            # check if they are already present:
            field = None
            new_atom_dict = {}
            with open(new_atomtypes) as file:
                for line in file:
                    if line.strip().startswith("["):
                        # Remove space and [ ], remove also comments
                        field = line.replace(" ", "").split("]")[0][1:]
                        continue

                    if (line[0] != ";" and line[0] != "#" and
                            line.strip() != ""):
                        # Remove commentary in the line
                        line_comment = line.split(';')
                        line_list = line_comment[0].split()

                        if field == 'atomtypes':
                            name, bond_type, mass, charge, ptype, sigma, \
                                epsilon = line_list[:7]
                            local_dict = {'bond_type': bond_type,
                                          'mass': float(mass),
                                          'charge': float(charge),
                                          'ptype': ptype,
                                          'sigma': float(sigma),
                                          'epsilon': float(epsilon)}
                            if name not in name_list:
                                new_atom_dict[name] = local_dict
                            else:
                                # Check if the new values are the same:
                                if local_dict != atom_dict[name]:
                                    logger.warning(
                                        'Atom types parameters for {} are '
                                        'different in {} and {}. Only one '
                                        'version is kept !!!'.format(
                                            name, atomtype_path,
                                            new_atomtypes))

            # Finally append the new param in the old one
            with open(atomtype_path, 'a') as file:
                for name, atom_dict in new_atom_dict.items():
                    file.write(' {:3}      {:3}         {:.5f}  {:.5f}'
                               '   {}     {:.5e}   {:.5e}\n'.format(
                                name,
                                atom_dict['bond_type'],
                                atom_dict['mass'],
                                atom_dict['charge'],
                                atom_dict['ptype'],
                                atom_dict['sigma'],
                                atom_dict['epsilon']))

    def add_intermolecular_restr(self, bond_list=[],
                                 angle_list=[], dihed_list=[]):
        """ Add inter molecular restraints in topologie file
        """

        self.inter_interact = True
        for bond in bond_list:
            ai, aj, funct, rA, kA, rB, kB = bond
            self.inter_bond_list.append({'ai': ai, 'aj': aj,
                                         'type': funct,
                                         'rA': rA, 'kA': kA,
                                         'rB': rB, 'kB': kB})
        for angle in angle_list:
            ai, aj, ak, funct, thA, kA, thB, kB = angle
            self.inter_angl_list.append({'ai': ai, 'aj': aj,
                                         'ak': ak, 'type': funct,
                                         'thA': thA, 'kA': kA,
                                         'thB': thB, 'kB': kB})
        for dihed in dihed_list:
            ai, aj, ak, al, funct, thA, kA, thB, kB = dihed
            self.inter_dihe_list.append({'ai': ai, 'aj': aj,
                                         'ak': ak, 'al': al,
                                         'type': funct,
                                         'thA': thA, 'kA': kA,
                                         'thB': thB, 'kB': kB})

    def change_mol_num(self, mol_name, mol_num):
        """ Update molecule number.
        And remove multiple molecule definition if they are consecutive.
        """

        for i, mol in enumerate(self.mol_comp):

            if mol['name'] == mol_name:
                mol['num'] = str(mol_num)
                # Find other mol_name and remove them:
                to_remove_index_list = []
                for j in range(i + 1, len(self.mol_comp)):
                    if self.mol_comp[j]['name'] == mol_name:
                        to_remove_index_list.append(j)
                    else:
                        break
                if len(to_remove_index_list) > 0:
                    for index in sorted(to_remove_index_list, reverse=True):
                        del self.mol_comp[index]
                break
        logger.info(self.mol_comp)

    def get_include_file_list(self):
        file_list = []
        # print("\n\nIncluded files: ")
        for itp in self.itp_list:
            if self.forcefield:
                if (self.forcefield['fullname'].split("/")[0] !=
                        itp.fullname.split("/")[0]):
                    # print(itp.path)
                    file_list.append(itp.path)
                    file_list = file_list + itp.get_include_file_list()
            else:
                file_list.append(itp.path)
                file_list = file_list + itp.get_include_file_list()

        return file_list

    def get_include_no_posre_file_list(self):
        file_list = []
        # print("\n\nIncluded files: ")
        for itp in self.itp_list:
            if self.forcefield:
                if (self.forcefield['fullname'].split("/")[0] !=
                        itp.fullname.split("/")[0]):
                    file_list.append(itp.path)
            else:
                file_list.append(itp.path)
        return file_list

    def copy_top_and_dependancies(self, dest_file):

        dest_folder = os_command.get_directory(dest_file)
        logger.info("Copy topologie file and dependancies")

        if self.path != os.path.abspath(dest_file):
            shutil_copy(self.path, os.path.abspath(dest_file))

        self.copy_dependancies(dest_folder)

    def copy_dependancies(self, dest_folder):

        dest_folder = os.path.abspath(dest_folder)

        file_to_copy = self.get_include_file_list()
        # print(file_to_copy)
        for file in file_to_copy:
            if os.path.abspath(os_command.get_directory(file)) !=\
                    os.path.abspath(dest_folder):
                # print("Copy: "+file+":to: "+dest_folder)
                shutil_copy(file, dest_folder)

    def change_mol_name(self, old_name, new_name):

        # Change name in mol_comp
        name_change_flag = False
        for mol in self.mol_comp:
            if mol['name'] == old_name:
                name_change_flag = True
                mol['name'] = new_name
        # Change in itp :
        for itp in self.itp_list:
            itp.change_mol_name(old_name, new_name)

        if name_change_flag:
            self.write_file(self.path)
