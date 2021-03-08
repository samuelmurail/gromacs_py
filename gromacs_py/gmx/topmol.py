#!/usr/bin/env python3
# coding: utf-8
##################################
# #######   GROMACS 5   ##########
##################################


import os
import logging
import sys

from os_command_py import os_command

from .rtp import Rtp

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


def show_debug():
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.DEBUG)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))


class TopMol:
    """Individual molecule topologie
    """

    def __init__(self, name, nrexcl):
        """Read an itp file and extract [atoms] field
        """
        self.name = name
        self.nrexcl = nrexcl
        self.atom_dict = {}
        self.bond_list = []
        self.cons_list = []
        self.pair_list = []
        self.angl_list = []
        self.dihe_list = []
        self.cmap_list = []
        self.vs3_list = []
        self.vs4_list = []
        self.pos_restr = []
        self.if_pos_restr = []

    def get_charge(self):
        self.charge = 0
        for atom in self.atom_dict.values():
            self.charge += atom['charge']
        return round(self.charge, 4)

    def get_res_num(self):
        self.res_num = -1
        for atom in self.atom_dict.values():
            if atom["res_num"] != self.res_num:
                self.res_num += 1
        return self.res_num

    def get_selection_index(self, selec_dict={'atom_name': ['CA']}):
        """Return the atom index to add posre
        """

        index_list = []
        for index, atom in self.atom_dict.items():
            selected = True
            for field, selection in selec_dict.items():
                if atom[field] not in selection:
                    selected = False
                    break
            if selected:
                index_list.append(index)

        return index_list

    def write_file(self, filout):
        filout.write("\n[ moleculetype ]\n; Name            nrexcl\n")
        filout.write("{}\t\t{}\n".format(self.name, self.nrexcl))
        # Print Atoms field
        filout.write("\n[ atoms ]\n;   nr       type  resnr residue  atom   "
                     "cgnr     charge       mass  typeB    chargeB      "
                     "massB\n")
        tot_charge = 0
        for atom in self.atom_dict.values():
            tot_charge += atom['charge']
            # print('ATOM:',atom)
            filout.write("{:>6}{:>11}{:>7}{:>7}{:>7}{:>7}{:>11.5f}{:>11} "
                         "  ; qtot {:<6.2f} \n".
                         format(atom['num'],
                                atom['atom_type'],
                                atom['res_num'],
                                atom['res_name'],
                                atom['atom_name'], atom['charge_num'],
                                atom['charge'],
                                atom['mass'],
                                tot_charge))
        # Print bonds field
        if self.bond_list:
            filout.write("\n[ bonds ]\n;  ai    aj funct            c0        "
                         "    c1            c2            c3\n")
            for param in self.bond_list:
                filout.write("{:>6}{:>6}{:>6}{:>13}{:>13}\n".format(
                    param['ai'], param['aj'], param['funct'],
                    param['r'], param['k']))
        # Print constraints field
        if self.cons_list:
            filout.write("\n[ constraints ]\n;  ai    aj funct            c0  "
                         "          c1\n")
            for param in self.cons_list:
                filout.write("{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['aj'], param['funct']))
        # Print pairs field
        if self.pair_list:
            filout.write("\n[ pairs ]\n;  ai    aj funct            c0        "
                         "    c1            c2            c3\n")
            for param in self.pair_list:
                filout.write("{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['aj'], param['funct']))
        # Print angles field
        if self.angl_list:
            filout.write("\n[ angles ]\n;  ai    aj    ak funct            c0 "
                         "           c1            c2            c3\n")
            for param in self.angl_list:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}\n".format(
                    param['ai'], param['aj'],
                    param['ak'], param['funct'],
                    param['theta'], param['cth']))
        # Print dihedrals field
        if self.dihe_list:
            filout.write("\n[ dihedrals ]\n;  ai    aj    ak    al funct      "
                         "      c0            c1            c2            c3  "
                         "          c4            c5\n")
            for param in self.dihe_list:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}"
                             "{:>6}\n".format(param['ai'], param['aj'],
                                              param['ak'], param['al'],
                                              param['funct'], param['phase'],
                                              param['kd'], param['pn']))
        # Print virtual_sites3 field
        if self.cmap_list:
            filout.write("\n[ cmap ]\n;  ai    aj    ak    al    am funct\n")
            for param in self.cmap_list:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['aj'],
                    param['ak'], param['al'],
                    param['am'], param['funct']))
        # Print virtual_sites3 field
        if self.vs3_list:
            filout.write("\n[ virtual_sites3 ]\n;  ai    aj    ak    al funct"
                         "            c0            c1\n")
            for param in self.vs3_list:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['aj'],
                    param['ak'], param['al'],
                    param['funct']))
        # Print virtual_sites3 field
        if self.vs4_list:
            filout.write("\n[ virtual_sites4 ]\n;  ai    aj    ak    al    am"
                         " funct            c0            c1            c2\n")
            for param in self.vs4_list:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['aj'],
                    param['ak'], param['al'],
                    param['am'], param['funct']))
        # Print position restraints
        if self.pos_restr:
            filout.write("\n[ position_restraints ]\n;  ai    funct         "
                         "kx          ky          kz\n")
            for param in self.pos_restr:
                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['funct'],
                    param['kx'], param['ky'],
                    param['kz']))
        if len(self.if_pos_restr) > 0:
            def_str = ''
            def_end_flag = False
            for param in self.if_pos_restr:

                if param['def'] != def_str:
                    def_str = param['def']
                    if def_end_flag:
                        filout.write("#endif\n\n")

                    def_end_flag = True
                    filout.write("\n#ifdef " + param['def'] + "\n")
                    filout.write("\n[ position_restraints ]\n;  ai    funct "
                                 "        kx          ky          kz\n")

                filout.write("{:>6}{:>6}{:>6}{:>6}{:>6}\n".format(
                    param['ai'], param['funct'],
                    param['kx'], param['ky'],
                    param['kz']))
            filout.write("#endif\n\n")

    def delete_atom(self, index_list):
        # Remove atom:
        for i in index_list:
            del self.atom_dict[i]

        # Create the dict to have all atom num consecutive staring from 1
        dict_atom_index = {}
        for i, atom in sorted(enumerate(self.atom_dict.items())):
            # print(i, atom)
            dict_atom_index[atom[0]] = i + 1

        # Modify all atom to have correct index
        new_atom_dict = {}
        for i, atom in self.atom_dict.items():
            local_atom = atom
            new_index = dict_atom_index[local_atom['num']]
            local_atom['num'] = new_index
            local_atom['charge_num'] = new_index
            new_atom_dict[new_index] = local_atom

        self.atom_dict = new_atom_dict

        #  Modify all bond, cons, pair ... to have correct index:
        new_bond_list = []
        for i, param in enumerate(self.bond_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']]})
                new_bond_list.append(param)
        self.bond_list = new_bond_list

        new_cons_list = []
        for i, param in enumerate(self.cons_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']]})
                new_cons_list.append(param)
        self.cons_list = new_cons_list

        new_pair_list = []
        for i, param in enumerate(self.pair_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']]})
                new_pair_list.append(param)
        self.pair_list = new_pair_list

        new_angl_list = []
        for i, param in enumerate(self.angl_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list) or
                    (param['ak'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']],
                              'ak': dict_atom_index[param['ak']]})
                new_angl_list.append(param)

        self.angl_list = new_angl_list

        new_dihe_list = []
        for i, param in enumerate(self.dihe_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list) or
                    (param['ak'] in index_list) or
                    (param['al'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']],
                              'ak': dict_atom_index[param['ak']],
                              'al': dict_atom_index[param['al']]})
                new_dihe_list.append(param)

        self.dihe_list = new_dihe_list

        new_vs3_list = []
        for i, param in enumerate(self.vs3_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list) or
                    (param['ak'] in index_list) or
                    (param['al'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']],
                              'ak': dict_atom_index[param['ak']],
                              'al': dict_atom_index[param['al']]})
                new_vs3_list.append(param)
        self.vs3_list = new_vs3_list

        new_cmap_list = []
        for i, param in enumerate(self.cmap_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list) or
                    (param['ak'] in index_list) or
                    (param['al'] in index_list) or
                    (param['am'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']],
                              'ak': dict_atom_index[param['ak']],
                              'al': dict_atom_index[param['al']],
                              'am': dict_atom_index[param['am']]})
                new_cmap_list.append(param)
        self.cmap_list = new_cmap_list

        new_vs4_list = []
        for i, param in enumerate(self.vs4_list):
            if not ((param['ai'] in index_list) or
                    (param['aj'] in index_list) or
                    (param['ak'] in index_list) or
                    (param['al'] in index_list) or
                    (param['am'] in index_list)):
                param.update({'ai': dict_atom_index[param['ai']],
                              'aj': dict_atom_index[param['aj']],
                              'ak': dict_atom_index[param['ak']],
                              'al': dict_atom_index[param['al']],
                              'am': dict_atom_index[param['am']]})
                new_vs4_list.append(param)
        self.vs4_list = new_vs4_list

    def correct_charge_type(self, forcefield, index_list=None):
        """ Correct the charge and atom type of an itp object,
        base on a ff `.rtp` file.
        This is specially usefull, to correct charge of termini resiudes
        of a cyclicized peptide.

        if index_list is None, will correct all atom charges, if not, only
        atoms listed in index_list.
        """

        # First extract charge and type from the ff .rtp file
        forcefield_path = forcefield['path']
        rtp_path = os.path.abspath(os.path.join(
            os_command.get_directory(forcefield_path), 'aminoacids.rtp'))
        if not os_command.check_file_exist(rtp_path):
            rtp_path = os.path.abspath(os.path.join(
                os_command.get_directory(forcefield_path), 'merged.rtp'))

        logger.info('Read rtp file : {}'.format(rtp_path))
        ff_rtp = Rtp(rtp_path)

        # Correct charges and type:
        # for atom_num, atom in self.atom_dict.items():
        if index_list is None:
            index_list = self.atom_dict.keys()

        for atom_num in index_list:

            atom = self.atom_dict[atom_num]

            res_name = atom['res_name']
            resid = atom['res_num']
            atom_name = atom['atom_name']

            # With Amber disulfide cys are called CYX in the rtp and CYS in
            # the itp.
            # Don't get sense but need to check the protonation state of cys
            # By counting atoms
            # With charmm disuldie cys is CYS2
            if res_name == 'CYS':
                if len(self.get_selection_index(
                        selec_dict={'res_num': [resid],
                                    'res_name': ['CYS']})) == 10:
                    if forcefield['name'].startswith('amber'):
                        res_name = 'CYX'
                    elif forcefield['name'].startswith('charmm'):
                        res_name = 'CYS2'
            # print(res_name)
            # With gromacs all histidine are called HIS !
            if res_name == 'HIS':
                if len(self.get_selection_index(
                        selec_dict={'res_num': [resid],
                                    'res_name': ['HIS']})) == 18:
                    if forcefield['name'].startswith('amber'):
                        res_name = 'HIP'
                    elif forcefield['name'].startswith('charmm'):
                        res_name = 'HSP'
                elif len(self.get_selection_index(
                        selec_dict={'res_num': [resid],
                                    'atom_name': ['HE2']})) == 1:
                    if forcefield['name'].startswith('amber'):
                        res_name = 'HIE'
                    elif forcefield['name'].startswith('charmm'):
                        res_name = 'HSE'
                else:
                    if forcefield['name'].startswith('amber'):
                        res_name = 'HID'
                    elif forcefield['name'].startswith('charmm'):
                        res_name = 'HSD'

            # print(atom)
            # print(ff_rtp.res_dict[res_name]['atom'][atom_name])

            if atom_name in ff_rtp.res_dict[res_name]['atom']:
                atom_type = \
                    ff_rtp.res_dict[res_name]['atom'][atom_name]['type']
                atom_charge = \
                    ff_rtp.res_dict[res_name]['atom'][atom_name]['charge']
                # print('Correct residue {:4} atom {:4} atom type {:4} '
                #       'to {:4}'.format(res_name, atom['atom_name'],
                #                        atom['atom_type'], atom_type))
                # print('Correct charge {:4} '
                #       'to {:4}'.format(self.atom_dict[atom_num]['charge'],
                #                        atom_charge))
                if atom_type != atom['atom_type']:
                    logger.warning('Correct residue {:4} atom {:4} atom '
                                   'type {:4} to {:4}'.format(
                                        res_name, atom['atom_name'],
                                        atom['atom_type'], atom_type))
            else:
                logger.warning('Can\'t find residue {:4} atom {:4} atom '
                               'type {:4} parameters in forcefield'.format(
                                    res_name, atom['atom_name'],
                                    atom['atom_type']))
            self.atom_dict[atom_num]['atom_type'] = atom_type
            self.atom_dict[atom_num]['charge'] = atom_charge
