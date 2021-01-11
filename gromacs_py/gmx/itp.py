#!/usr/bin/env python3
# coding: utf-8
##################################
# #######   ITP Class   ##########
##################################


import os
import logging
import sys

from os_command_py import os_command

from .topmol import TopMol

# Logging
logger = logging.getLogger(__name__)


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
    from . import topmol
    topmol.show_log()


# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "1.2.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Prototype"


# Position restraints files:
def write_index_posre_file(atom_index_list, posre_file, type_val=1,
                           fc=[1000, 1000, 1000]):
    """Write a pos restraint file based on atom index list
    """

    filout = open(posre_file, 'w')
    filout.write("; Position restraint file created by " + __name__ + "\n")
    filout.write("[ position_restraints ]\n")
    filout.write("; atom  type      fx      fy      fz\n")

    for index in atom_index_list:
        filout.write("{:6d}{:6d}{:6d}{:6d}{:6d} \n".format(
            index, type_val, fc[0], fc[1], fc[2]))

    filout.write("\n")
    filout.close()


class Itp:
    """Itp topologie in gromacs format
    May contain several molecule definition, so
    itp class is a collection of top_mol object
    which are indivudial molecule topologies
    """

    def __init__(self, name, fullname, path):
        """Read an itp file and extract all the top_mol
        """

        # name is useless to fix, maybe only put shortpath and fullpath
        self.name = name
        self.fullname = fullname
        self.path = path
        self.top_mol_list = []
        self.posres_file = []
        self.read_file()

    def read_file(self):

        field = None
        ifdef = False
        posre_def = ""
        local_top = None
        local_top_flag = False

        with open(self.path) as file:
            for line in file:
                # print("Itp line: \"{}\" ".format(line))
                # Check posres include:
                if line[:6] == '#ifdef':
                    ifdef = True
                    line_list = line.split()
                    posre_def = line_list[1]
                if line[:6] == '#endif':
                    ifdef = False
                if ifdef and line[:8] == '#include':
                    line_list = line.split()
                    posre_file = line_list[1][1:-1]
                    self.posres_file.append(
                        {'def': posre_def, 'file': posre_file})
                if line.strip().startswith("["):
                    # Remove space and [ ], remove also comments
                    field = line.replace(" ", "").split("]")[0][1:]
                    continue
                if line[0] != ";" and line[0] != "#" and line.strip() != "":
                    # Remove commentary in the line
                    line_comment = line.split(';')
                    line_list = line_comment[0].split()
                    # print(line_list)

                    if field == 'moleculetype':
                        # Check if a top_mol already exist, if yes append it
                        # to the top_mol_list
                        if local_top_flag:
                            # print("Add mol topologie",local_top.name)
                            self.top_mol_list.append(local_top)
                        # Create a new top_mol
                        local_top = TopMol(line_list[0], line_list[1])
                        local_top_flag = True

                    elif field == 'atoms':
                        atom_num = int(line_list[0])
                        atom_type = line_list[1]
                        # In rare cases the residue num is not an integer
                        # With duplicate residues like : eg "184A"
                        try:
                            res_num = int(line_list[2])
                        except ValueError:
                            res_num = line_list[2]
                        res_name = line_list[3]
                        atom_name = line_list[4]
                        charge_num = int(line_list[5])
                        charge = float(line_list[6])
                        if len(line_list) > 7:
                            mass = float(line_list[7])
                        else:
                            mass = None
                        atom = {"num": atom_num,
                                "atom_type": atom_type,
                                "atom_name": atom_name,
                                "res_num": res_num,
                                "res_name": res_name,
                                "charge_num": charge_num,
                                "charge": charge,
                                "mass": mass}
                        local_top.atom_dict[atom_num] = atom
                    # Read only indexes not the parameters (c0, c1 ...)
                    # May need some modifications :
                    elif field == 'bonds':
                        # for i in range(3):
                        #    print("val:{} #{}#".format(i, line[i*5:(i+1)*5]))
                        # ai, aj, funct = [int(line[i*5:(i+1)*5]) for i in
                        # range(3)]
                        ai, aj, funct = [int(col) for col in line_list[:3]]
                        if len(line_list) < 5:
                            r = k = ''
                        else:
                            r, k = line_list[3:5]
                        local_top.bond_list.append({'ai': ai, 'aj': aj,
                                                    'funct': funct, 'r': r,
                                                    'k': k})
                    elif field == 'constraints':
                        ai, aj, funct = [int(col) for col in line_list[:3]]
                        local_top.cons_list.append({'ai': ai, 'aj': aj,
                                                    'funct': funct})
                    elif field == 'pairs':
                        ai, aj, funct = [int(col) for col in line_list[:3]]
                        local_top.pair_list.append({'ai': ai, 'aj': aj,
                                                    'funct': funct})
                    elif field == 'angles':
                        ai, aj, ak, funct = [int(col) for col in line_list[:4]]
                        if len(line_list) < 6:
                            theta = cth = ''
                        else:
                            theta, cth = line_list[4:6]
                        local_top.angl_list.append({'ai': ai, 'aj': aj,
                                                    'ak': ak, 'funct': funct,
                                                    'theta': theta,
                                                    'cth': cth})
                    elif field == 'dihedrals':
                        ai, aj, ak, al, funct = [int(col) for col in
                                                 line_list[:5]]
                        if len(line_list) < 7:
                            phase = kd = pn = ''
                        else:
                            phase, kd = line_list[5:7]
                            pn = int(line_list[7])
                        local_top.dihe_list.append({'ai': ai, 'aj': aj,
                                                    'ak': ak, 'al': al,
                                                    'funct': funct,
                                                    'phase': phase, 'kd': kd,
                                                    'pn': pn})
                    elif field == 'virtual_sites3':
                        ai, aj, ak, al, funct = [
                            int(col) for col in line_list[:5]]
                        local_top.vs3_list.append({'ai': ai, 'aj': aj,
                                                   'ak': ak, 'al': al,
                                                   'funct': funct})
                    elif field == 'cmap':
                        ai, aj, ak, al, am, funct = [
                            int(col) for col in line_list[:6]]
                        local_top.cmap_list.append({'ai': ai, 'aj': aj,
                                                    'ak': ak, 'al': al,
                                                    'am': am, 'funct': funct})
                    elif field == 'virtual_sites4':
                        ai, aj, ak, al, am, funct = [
                            int(col) for col in line_list[:6]]
                        local_top.vs4_list.append({'ai': ai, 'aj': aj,
                                                   'ak': ak, 'al': al,
                                                   'am': am, 'funct': funct})
                    elif field == 'position_restraints':
                        # With MARTINI topologie kx, ky, kz can be strings,
                        # so no int() conversions
                        ai, funct, kx, ky, kz = [col for col in line_list[:5]]
                        if ifdef:
                            local_top.if_pos_restr.append(
                                {'def': posre_def, 'ai': int(ai),
                                 'funct': int(funct), 'kx': kx,
                                 'ky': ky, 'kz': kz})
                        else:
                            local_top.pos_restr.append(
                                {'ai': int(ai), 'funct': int(funct), 'kx': kx,
                                 'ky': ky, 'kz': kz})

                    # else:
                    #   raise ValueError('Unknown field : '+field)

        # Needed for empty topologies like aditional ff parameters:
        if local_top_flag:
            self.top_mol_list.append(local_top)

    def write_file(self, itp_file):
        filout = open(itp_file, 'w')
        filout.write("; Itp file created by " + __name__ + "\n")

        for top_mol in self.top_mol_list:
            logger.info(top_mol.name)
            top_mol.write_file(filout)
        for posre in self.posres_file:
            filout.write("\n#ifdef " + posre['def'] + "\n")
            filout.write("#include \"" + posre['file'] + "\"\n")
            filout.write("#endif\n\n")

        filout.close()

    def display(self):
        logger.info('- ITP file: {}'.format(self.name))
        logger.info("- molecules defined in the itp file:")
        for top_mol in self.top_mol_list:
            logger.info("* {}".format(top_mol.name))

    def charge(self, mol_name):
        for top_mol in self.top_mol_list:
            # print(mol_name, top_mol.name)
            if top_mol.name == mol_name:
                return top_mol.get_charge()
        return None

    def res_num(self, mol_name):
        for top_mol in self.top_mol_list:
            # print(mol_name, top_mol.name)
            if top_mol.name == mol_name:
                return top_mol.get_res_num()
        return None

    def add_posre(self, mol_name, posre_name, selec_dict, fc, replace=True):
        for top_mol in self.top_mol_list:
            # print(mol_name, top_mol.name)
            if top_mol.name == mol_name:
                index_posre = top_mol.get_selection_index(
                    selec_dict=selec_dict)
                if index_posre:
                    # Create the posre itp file :
                    logger.debug("Posre for : {}".format(top_mol.name))
                    posre_file_name = os.path.abspath(
                        os.path.join(
                            os_command.get_directory(self.path),
                            self.name + "_posre_" + posre_name + ".itp"))
                    # posre_file_name = self.name+"_posre_"+posre_name+".itp"
                    write_index_posre_file(atom_index_list=index_posre,
                                           posre_file=posre_file_name,
                                           type_val=1, fc=fc)
                    # Add the posre include in the mol itp file:
                    # Need to solve the problem of #ifdef location with .top
                    # files containing self top
                    if replace:
                        # Remove previous def:
                        new_content = ""
                        with open(self.path) as file:
                            posre_def = None
                            for line in file:
                                # print("Itp line: \"{}\" ".format(line))
                                # Check posres include:
                                if line[:6] == '#ifdef':
                                    line_list = line.split()
                                    posre_def = line_list[1]

                                if posre_def != posre_name:
                                    new_content += line

                                if line[:6] == '#endif':
                                    posre_def = None
                        # Add the new one:
                        new_content += '#ifdef ' + posre_name + '\n'
                        new_content += '#include \"' +\
                            os.path.basename(posre_file_name) + "\" \n"
                        new_content += '#endif \n\n'
                        # Save the file:
                        file = open(self.path, "w")
                        file.write(new_content)
                        file.close()

                    else:
                        with open(self.path, 'a') as file:
                            file.write('#ifdef ' + posre_name + '\n')
                            # file.write('#include \"'+self.name+"_posre_"+\
                            # posre_name+".itp\" \n")
                            file.write(
                                '#include \"' +
                                os.path.basename(posre_file_name) + "\" \n")
                            file.write('#endif \n\n')

    def set_top_mol_name(self, new_name):
        if len(self.top_mol_list) == 1:
            self.top_mol_list[0].name = new_name
        else:
            raise ValueError('Cannot set top mol name with multiple top mol '
                             'in itp')

    def get_include_file_list(self):
        file_list = []
        for posre in self.posres_file:
            # print(posre)
            file_list.append(os.path.abspath(os.path.join(
                os_command.get_directory(self.path), posre['file'])))
        return file_list

    def change_mol_name(self, old_name, new_name):

        name_change_flag = False
        for top_mol in self.top_mol_list:
            if top_mol.name == old_name:
                name_change_flag = True
                top_mol.name = new_name
        if name_change_flag:
            self.write_file(self.path)
