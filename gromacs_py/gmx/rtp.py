#!/usr/bin/env python3
# coding: utf-8
##################################
# #######   RTP Class   ##########
##################################


class Rtp:
    """Individual molecule topologie
    """

    def __init__(self, path):
        """Read an itp file and extract [atoms] field
        """
        self.path = path
        self.res_dict = {}
        self.read_file()

        # Then assigne the correct atom charge and type

    def read_file(self):

        field = None
        atom_dict = {}
        bond_list = []
        impr_list = []
        res_name = None

        with open(self.path) as file:
            for line in file:
                line_strip = line.strip()
                if len(line_strip) > 0:
                    if line_strip[0] == "[":
                        # Remove space and [ ], remove also comments
                        field = line.replace(" ", "").split("]")[0][1:]

                        if field not in ['bondedtype', 'atoms', 'bonds',
                                         'impropers', '']:

                            if len(atom_dict) > 0:
                                residue = {}
                                residue['atom'] = atom_dict
                                residue['bond'] = bond_list
                                residue['impr'] = impr_list
                                self.res_dict[res_name] = residue

                            res_name = field
                            atom_dict = {}
                            bond_list = []
                            impr_list = []
                        continue

                    if line_strip[0] != ";" and line_strip[0] != "#":
                        line_list = line_strip.split()

                        if field == 'atoms':
                            atom_name = line_list[0]
                            atom_type = line_list[1]
                            atom_charge = float(line_list[2])
                            charge_group = int(line_list[3])
                            atom_dict[atom_name] = {
                                "type": atom_type,
                                "charge": atom_charge,
                                "charge_group": charge_group}

                        elif field == 'bonds':
                            bond_list.append({'ai': line_list[0],
                                              'aj': line_list[1]})

                        elif field == 'impropers':
                            impr_list.append({'ai': line_list[0],
                                              'aj': line_list[1],
                                              'ak': line_list[2],
                                              'al': line_list[3]})
