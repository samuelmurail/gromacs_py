#!/usr/bin/env python3

""" Minimize a pdb file structure
"""

import argparse
import shutil
from gromacs_py import gmx
from pdb_manip_py import pdb_manip

__author__ = "Samuel Murail"


def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(
        description='Minimize a cyclic peptide structure in 2 steps, the '
        'first step without bonds constraints and the second step with '
        'bonds constraints')
    parser.add_argument('-f', action="store", dest="f",
                        help='Input PDB file', type=str, required=True)
    parser.add_argument('-n', action="store", dest="name",
                        help='Output file name', type=str, required=True)
    parser.add_argument('-dir', action="store", dest="out_dir",
                        help='Output directory for intermediate files',
                        type=str, default="./tmp_em")
    parser.add_argument('-m_steps', action="store", dest="min_steps",
                        help='Minimisation nsteps, default=1000', type=int,
                        default=1000)
    parser.add_argument('-keep', action="store_true", dest="keep_flag",
                        help='Flag to keep temporary files (without flag '
                        'output directory will be delete')
    parser.add_argument('-cyclic', action="store_true", dest="cyclic_flag",
                        help='Flag to indicate if the peptide/protein is '
                        'cyclic')
    parser.add_argument('-nt', action="store", dest="nt",
                        help='Total number of threads to start, default=0',
                        type=float, default=0)
    parser.add_argument('-ntmpi', action="store", dest="ntmpi",
                        help='Number of thread-MPI threads to start, '
                        'default=0', type=float, default=0)
    parser.add_argument('-gpu_id', action="store", dest="gpuid",
                        help='List of GPU device id-s to use, default=\"\" ',
                        default="None")
    parser.add_argument('-keep_segid', action="store_true", dest="keep_segid",
                        help='Flag to indicate if the original chain/segid '
                        'should be kept')
    parser.add_argument('-add_ter', action="store_true", dest="add_termini",
                        help='Flag to indicate if TER line should be included'
                        ' between chains and if residues in the input pdb '
                        'file have non-consecutive residues')

    return parser


def get_structure_string(pdb_file):
    """Return a pdb_file as a pdb string.
    Add TER line between chains or between
    non consecutive residues.
    """

    str_out = ""
    ref_coor = pdb_manip.Coor(pdb_file)

    _, first_atom = list(sorted(ref_coor.atom_dict.items()))[0]
    chain = first_atom['chain']
    res_num = first_atom['res_num']

    for atom_num, atom in sorted(ref_coor.atom_dict.items()):
        # Atom name should start a column 14, with the type of atom ex:
        #   - with atom type 'C': ' CH3'
        # for 2 letters atom type, it should start at coulumn 13 ex:
        #   - with atom type 'FE': 'FE1'
        name = atom["name"]
        if len(name) <= 3 and name[0] in ['C', 'H', 'O', 'N', 'S', 'P']:
            name = " " + name

        # Add the TER
        if chain != atom["chain"]:
            print('add TER for chain {}'.format(atom['chain']))
            str_out += "TER\n"
            chain = atom["chain"]
            res_num = atom['res_num']
        elif res_num != atom['res_num']:
            if res_num == atom['res_num']-1:
                res_num += 1
            else:
                print('add TER between res {} and {} of chain {}'.format(
                    res_num, atom['res_num'], atom['chain']))
                res_num = atom['res_num']
                str_out += "TER\n"

        # Note : Here we use 4 letter residue name.
        str_out += "{:6s}{:5d} {:4s}{:1s}{:4s}{:1s}{:4d}{:1s}"\
                   "   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"\
                   "          {:2s}\n".format(
                        atom["field"],
                        atom["num"],
                        name,
                        atom["alter_loc"],
                        atom["res_name"],
                        atom["chain"],
                        atom["res_num"],
                        atom["insert_res"],
                        atom["xyz"][0],
                        atom["xyz"][1],
                        atom["xyz"][2],
                        atom["occ"],
                        atom["beta"],
                        atom['elem'])

    return str_out


def add_ter(pdb_file, pdb_out):
    """ Add TER between chain or
    when residue number are not consecutive.
    """

    coor_str = get_structure_string(pdb_file)

    filout = open(pdb_out, 'w')
    filout.write(coor_str)
    filout.close()


def get_chain_res_list(pdb_file):
    """ Get chain and residue list
    In order to remap them to gromacs
    pdb.
    """

    ref_coor = pdb_manip.Coor(pdb_file)
    ca_coor = ref_coor.select_part_dict(selec_dict={'name': ['CA']})

    chain_res_list = []

    for atom_num, atom in sorted(ca_coor.atom_dict.items()):
        chain_res_list.append([atom['chain'], atom['res_num']])

    return(chain_res_list)


def set_chain_res_list(pdb_file, chain_res_list):
    """ Set the chain and residue as in the provided list.
    Save the new pdb and overwrite the pdb input.
    """

    ref_coor = pdb_manip.Coor(pdb_file)
    i = 0

    for atom_num, atom in sorted(ref_coor.atom_dict.items()):

        if atom['res_num'] == chain_res_list[i][1]:
            atom['chain'] = chain_res_list[i][0]
        elif (i + 1 >= len(chain_res_list)):
            print('Mismatch probably, ligand is included:',
                  chain_res_list[i], atom)
            break
        elif atom['res_num'] == chain_res_list[i + 1][1]:
            i += 1
            atom['chain'] = chain_res_list[i][0]
        else:
            print('Mismatch :', chain_res_list[i],
                  chain_res_list[i + 1], atom)
            print('\n' * 10 + 'WRONG' + '\n' * 10)
            break

    ref_coor.write_pdb(pdb_file, check_file_out=False)


if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()

    if args.add_termini:
        add_ter(args.f, args.f[:-4] + '_good_ter.pdb')
        input_pdb = args.f[:-4] + '_good_ter.pdb'
    else:
        input_pdb = args.f

    vsite = "none"
    peptide = gmx.GmxSys(name=args.name, coor_file=input_pdb)
    peptide.nt = args.nt
    peptide.ntmpi = args.ntmpi
    if args.gpuid != "None":
        peptide.gpu_id = args.gpuid

    if args.cyclic_flag:
        peptide.cyclic_peptide_top(out_folder=args.out_dir + '/top')
    else:
        peptide.add_top(out_folder=args.out_dir + '/top',
                        pdb2gmx_option_dict={'vsite': 'no',
                                             'ignh': 'yes',
                                             'ter': 'no'})

    peptide.em(out_folder=args.out_dir + '/em', name="min_" + args.name,
               nsteps=args.min_steps, posres="", nstxout=1000,
               create_box_flag=True,
               constraints="none")

    peptide.convert_trj(traj=False)

    # Get the minimised structure:
    shutil.copyfile(peptide.coor_file, args.name + '.pdb')

    if args.keep_segid:
        chain_res_list = get_chain_res_list(input_pdb)
        set_chain_res_list(args.name + '.pdb', chain_res_list)

    # Keep or not the intermediate files:
    if not args.keep_flag:
        shutil.rmtree(args.out_dir, ignore_errors=True)
    else:
        peptide.display()
