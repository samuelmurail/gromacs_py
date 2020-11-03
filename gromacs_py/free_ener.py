#!/usr/bin/env python3
# coding: utf-8
##################################
# ######   Free Energy   #########
##################################

import logging
import math
import os
import copy
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from os_command_py import os_command
from pdb_manip_py import pdb_manip

from . import gmx
from .gmx import TopSys, GmxSys, GROMACS_MOD_DIRNAME
from .tools import monitor


def show_log(pdb_manip_log=True):
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.INFO)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))

    gmx.show_log()
    if pdb_manip_log:
        # Show pdb_manip Logs:
        pdb_manip.show_log()


# Logging
logger = logging.getLogger(__name__)


# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Development"


# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    logger.info("Gromacs cannot be found")
    GMX_BIN = ""
else:
    GMX_BIN = os_command.which('gmx')


# Boltzmann Constant
KB = 8.31446261815324


class FreeEner:

    def __init__(self, mol_name, out_folder, unit='kcal'):
        self.mol_name = mol_name
        self.unit = unit
        self.out_folder = out_folder
        # Lambda
        self.lambda_coul = []
        self.lambda_vdw = []
        self.lambda_restr = []
        # Xvg file list
        self.xvg_file_list = []
        self.temp = None
        self.smile = None

    @property
    def conv_fac(self):
        """ Conversion factor from kT to
        self.unit
        """
        if self.unit == 'kcal':
            return(0.593)
        elif self.unit == 'kJ':
            return(4.184 * 0.593)
        elif self.unit == 'kT':
            return(1)
        elif self.unit == 'logP':
            return(-1.365679 * KB * self.temp / 1000)

    @property
    def unit_name(self):
        if self.unit == 'kcal':
            return('kcal/mol')
        elif self.unit == 'kJ':
            return('kJ/mol')
        elif self.unit == 'kT':
            return('kT')
        elif self.unit == 'logP':
            return('logP')

    @property
    def unit_graph(self):
        if self.unit == 'kcal':
            return(r'kcal \cdot mol^{-1}')
        elif self.unit == 'kJ':
            return(r'kJ \cdot mol^{-1}')
        elif self.unit == 'kT':
            return(r'k_{B}T')
        elif self.unit == 'logP':
            return(r'log \; P')

    @staticmethod
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

    def water_box_from_SMILE(self, smile, method_3d='rdkit',
                             iter_num=5000):
        """ Create water box with a molecule
        """

        self.smile = smile

        os.makedirs(self.out_folder, exist_ok=True)

        # Create 3d structure from SMILE
        pdb_file = os.path.join(self.out_folder, f'{self.mol_name}.pdb')
        charge = FreeEner.smile_to_pdb(smile, pdb_file,
                                       method_3d=method_3d,
                                       mol_name=self.mol_name,
                                       iter_num=iter_num)
        logger.info(f'ligand charge is {charge}')

        # Topologie and system creation
        self.gmxsys = GmxSys(name=self.mol_name, coor_file=pdb_file)
        self.gmxsys.prepare_top_ligand(
            out_folder=os.path.join(self.out_folder, 'mol_top'),
            ff='amber99sb-ildn', include_mol={self.mol_name: smile})

        self.gmxsys.create_box(name=None, dist=1.1)
        self.gmxsys.solvate_box(
            out_folder=os.path.join(self.out_folder, 'mol_water_top'))

    def equilibrate_solvent_box(self, em_steps=10000, dt=0.002, prod_time=10,
                                temp=300):
        """ Create water box with a molecule
        """
        # EM
        self.gmxsys.em_2_steps(out_folder=os.path.join(
                                self.out_folder, 'mol_water_em'),
                               no_constr_nsteps=em_steps,
                               constr_nsteps=em_steps,
                               posres="",
                               emtol=0.1, nstxout=100, maxwarn=3,
                               create_box_flag=False)
        mdp_options = {'tc-grps': 'System',
                       'tau_t': '0.1',
                       'ref_t': temp,
                       }
        # Equilibration
        # first short step with very small dt
        self.gmxsys.production(
            out_folder=os.path.join(self.out_folder, 'mol_water_prod'),
            name='short_equi',
            nsteps=50000, dt=0.0005,
            maxwarn=1, **mdp_options)
        # second step
        prod_steps = 1000 * prod_time / dt
        self.gmxsys.production(
            out_folder=os.path.join(self.out_folder, 'mol_water_prod'),
            nsteps=prod_steps,
            dt=dt, maxwarn=1, **mdp_options)

    def prepare_complex_pdb(self, pdb_in, smile, ff='amber99sb-ildn'):
        """
        """
        self.gmxsys = GmxSys(name=self.mol_name, coor_file=pdb_in)
        self.gmxsys.prepare_top(out_folder=os.path.join(
                                    self.out_folder, 'complex_top'),
                                ff=ff, include_mol={self.mol_name: smile})
        self.gmxsys.solvate_add_ions(out_folder=os.path.join(
                                        self.out_folder, 'sys_top'),
                                     ion_C=0.15, maxwarn=3,
                                     create_box_flag=True)

    def equilibrate_complex(self, em_steps=10000, HA_time=0.25,
                            CA_time=0.5, CA_LOW_time=1.0,
                            dt=0.002, dt_HA=0.001, temp=300):
        """
        """

        self.gmxsys.em_2_steps(out_folder=os.path.join(
                                self.out_folder, 'prot_em'),
                               no_constr_nsteps=em_steps,
                               constr_nsteps=em_steps,
                               posres="",
                               maxwarn=3,
                               create_box_flag=False)

        self.gmxsys.convert_trj(select='0' + '\n System', center='yes')

        self.ref_coor = self.gmxsys.coor_file

        # Short equilibration with molecule in its own temp group
        mdp_options = {'nsteps': 20000,
                       'define': '-DPOSRES', 'dt': 0.0005,
                       'tc-grps': 'Protein {} Water_and_ions'.format(
                            self.mol_name),
                       'tau_t': '0.1 0.1 0.1',
                       'ref_t': f'{temp} {temp} {temp}'}

        equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                         "template/equi.mdp")

        self.gmxsys.run_md_sim(out_folder=os.path.join(
                                self.out_folder, "sys_short_equi"),
                               name="equi_HA_short_" + self.mol_name,
                               mdp_template=equi_template_mdp,
                               mdp_options=mdp_options, maxwarn=3,
                               pdb_restr=self.ref_coor)

        self.gmxsys.equi_three_step(out_folder=os.path.join(
                                     self.out_folder, 'sys_equi'),
                                    nsteps_HA=1000 * HA_time / dt_HA,
                                    nsteps_CA=1000 * CA_time / dt,
                                    nsteps_CA_LOW=1000 * CA_LOW_time / dt,
                                    dt=dt, dt_HA=dt_HA, maxwarn=3,
                                    pdb_restr=self.ref_coor)

    def run(self, lambda_coul_list, lambda_vdw_list, lambda_restr_list=[],
            mbar=False, dir_name='free_ener_run',
            em_steps=5000, nvt_time=10, npt_time=10, prod_time=100,
            dt=0.002, temp=300.0,
            temp_groups='Protein non-Protein', maxwarn=1,
            monitor_tool=monitor.PROGRESS_BAR):
        """Compute free energy to transfer a molecule from the
        system to vacum.

        :param out_folder: path of the output folder
        :type out_folder: str

        :param mol_name: Name of the molecule
        :type mol_name: str

        :param lambda_coul_list: List lambda points for Coulomb
        :type lambda_coul_list: list

        :param lambda_vdw_list: List lambda points for Lennard Jones
        :type lambda_vdw_list: list

        :param lambda_bond_list: List lambda points for restraints
        :type lambda_bond_list: list, default=[]

        :param mbar: MBAR flag
        :type mbar: bool, default=False

        :param em_steps: number of minimisation steps
        :type em_steps: int, default=5000

        :param nvt_time: Time (ps) of NVT equilibration
        :type nvt_time: int, default=10 ps

        :param npt_time:  Time (ps) of NPT equilibration
        :type npt_time: int, default=10 ps

        :param prod_time: Time (ps) of production run
        :type prod_time: float, default=100 ps

        :param dt: integration time step
        :type dt: float, default=0.002

        :param name: name of the simulation to run
        :type name: str, default=None

        :param temp: Temperature K
        :type temp: float, default=300.0

        :param temp_groups: Group(s) for temperature coupling
        :type temp_groups: str, default='Protein non-Protein'

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type rerun: dict, default=None

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.tpr
            * self.sim_name
            * self.coor_file
            * self.xtc

        """

        if monitor.isnotebook():
            from tqdm.notebook import tqdm
        else:
            from tqdm import tqdm

        # Remove lambda=0 for vdw and elec
        # If previous restr lambdas are computed
        if len(lambda_restr_list) > 0:
            lambda_coul_list = list(lambda_coul_list)
            lambda_coul_list.remove(0.0)
        if len(lambda_restr_list) + len(lambda_coul_list) > 0:
            lambda_vdw_list = list(lambda_vdw_list)
            lambda_vdw_list.remove(0.0)

        lambda_restr_num = len(lambda_restr_list)
        lambda_coul_num = len(lambda_coul_list)
        lambda_vdw_num = len(lambda_vdw_list)

        restr_lambdas = "".join([
            '{:.2f} '.format(i) for i in lambda_restr_list]) +\
            "1.00 " * (lambda_vdw_num + lambda_coul_num)
        coul_lambdas = "0.00 " * lambda_restr_num + "".join([
            '{:.2f} '.format(i) for i in lambda_coul_list]) +\
            "1.00 " * (lambda_vdw_num)
        vdw_lambdas = "0.00 " * (lambda_coul_num + lambda_restr_num) +\
            "".join(['{:.2f} '.format(i) for i in lambda_vdw_list])

        logger.info(f"Coulomb lambda : {coul_lambdas}\n"
                    f"Vdw lambda     : {vdw_lambdas}\n"
                    f"restr_lambdas  : {restr_lambdas}")

        # https://events.prace-ri.eu/event/674/attachments/618/896/MD_FreeEnergyTutorial.pdf
        # Suggest to use tau_t = 1.0 to avoid
        # over-damping the dynamics of water
        tau_t_str = " ".join(['{}'.format(1.0)
                              for _ in range(len(temp_groups.split()))])
        ref_t_str = " ".join(['{}'.format(temp)
                              for _ in range(len(temp_groups.split()))])

        free_ener_option_md = {'integrator': 'sd',
                               'dt': dt,
                               'constraints': 'all-bonds',
                               'nstcalcenergy': 50,
                               'tcoupl': '',
                               'tc_grps': temp_groups,
                               'tau_t': tau_t_str,
                               'ref_t': ref_t_str,
                               'vdwtype': 'cut_off',
                               'vdw_modifier': 'force_switch',
                               'rvdw_switch': 1.0,
                               'rvdw': 1.1,
                               'coulombtype': 'pme',
                               'rcoulomb': 1.1,
                               'free_energy': 'yes',
                               'init_lambda-state': 0,
                               'calc-lambda-neighbors': 1,
                               'delta_lambda': 0,
                               'coul_lambdas': coul_lambdas,
                               'vdw_lambdas': vdw_lambdas,
                               'bonded_lambdas': restr_lambdas,
                               'sc_alpha': 0.5,
                               'sc_power': 1,
                               'sc_sigma': 0.3,
                               'couple_moltype': self.mol_name,
                               'couple_lambda0': 'vdw-q',
                               'couple_lambda1': 'none',
                               'couple_intramol': 'no',
                               'nstdhdl': 50,
                               'separate_dhdl_file': 'yes'}

        self.xvg_file_list = []

        nvt_steps = int(nvt_time / dt)
        npt_steps = int(npt_time / dt)
        prod_steps = int(prod_time / dt)

        tot_step = (lambda_restr_num + lambda_coul_num + lambda_vdw_num) * (
            em_steps + nvt_steps + npt_steps + prod_steps)
        pbar = tqdm(total=tot_step)

        for i in range(lambda_restr_num + lambda_coul_num + lambda_vdw_num):

            logger.info('Compute lambda {} / {}'.format(
                i, lambda_restr_num + lambda_coul_num + lambda_vdw_num))

            xvg_file = FreeEner.compute_lambda_point(self.gmxsys, i,
                                                     self.mol_name,
                                                     os.path.join(
                                                        self.out_folder,
                                                        dir_name),
                                                     free_ener_option_md, pbar,
                                                     mbar=mbar,
                                                     em_steps=em_steps,
                                                     nvt_steps=nvt_steps,
                                                     npt_steps=npt_steps,
                                                     prod_steps=prod_steps,
                                                     dt=dt,
                                                     maxwarn=maxwarn,
                                                     monitor_tool=monitor_tool)
            self.xvg_file_list.append(xvg_file)

        self.temp = temp
        self.lambda_coul = list(lambda_coul_list)
        self.lambda_vdw = list(lambda_vdw_list)
        self.lambda_restr = list(lambda_restr_list)
        self.prod_time = prod_time

    @staticmethod
    def compute_lambda_point(gmx_sys, lambda_iter, mol_name, out_folder,
                             free_ener_option, pbar, mbar,
                             em_steps, nvt_steps, npt_steps,
                             prod_steps, dt=0.002, maxwarn=1,
                             monitor_tool=monitor.PROGRESS_BAR):
        """
        """

        start_sys = copy.deepcopy(gmx_sys)

        mini_template_mdp = os.path.join(
            GROMACS_MOD_DIRNAME, "template/mini.mdp")
        equi_template_mdp = os.path.join(
            GROMACS_MOD_DIRNAME, "template/equi.mdp")

        free_ener_option_md = copy.deepcopy(free_ener_option)
        free_ener_option_em = copy.deepcopy(free_ener_option)

        # Remove useless options for EM
        del(free_ener_option_em['integrator'],
            free_ener_option_em['dt'],
            free_ener_option_em['tc_grps'],
            free_ener_option_em['tau_t'],
            free_ener_option_em['ref_t'])

        sys_name = '{}_vdwq_{:02d}'.format(mol_name, lambda_iter)

        free_ener_option_em.update({'init_lambda-state': lambda_iter,
                                    'nsteps': em_steps})

        # Mini 5000 steps
        mol_sys = copy.deepcopy(start_sys)
        mol_sys.run_md_sim(out_folder=os.path.join(
                            out_folder, '00_em'),
                           name='em_' + sys_name,
                           mdp_template=mini_template_mdp,
                           mdp_options=free_ener_option_em, maxwarn=maxwarn,
                           monitor_tool=monitor_tool)
        pbar.update(em_steps)

        # MD
        free_ener_option_md.update({'init_lambda-state': lambda_iter})

        # NVT equilibration
        free_ener_option_md.update({'nsteps': nvt_steps,
                                    'pcoupl': 'no'})
        mol_sys.run_md_sim(out_folder=os.path.join(out_folder,
                                                   '01_equi_nvt'),
                           name='nvt_' + sys_name,
                           mdp_template=equi_template_mdp,
                           mdp_options=free_ener_option_md, maxwarn=maxwarn,
                           monitor_tool=monitor_tool)
        pbar.update(nvt_steps)

        # NPT equilibration
        free_ener_option_md.update({'nsteps': npt_steps,
                                    'pcoupl': 'parrinello-Rahman'})
        mol_sys.run_md_sim(out_folder=os.path.join(out_folder,
                                                   '02_equi_npt'),
                           name='npt_' + sys_name,
                           mdp_template=equi_template_mdp,
                           mdp_options=free_ener_option_md, maxwarn=maxwarn,
                           monitor_tool=monitor_tool)
        pbar.update(npt_steps)

        # Production
        if mbar:
            free_ener_option_md.update({'calc-lambda-neighbors': -1})

        free_ener_option_md.update({'nsteps': prod_steps})
        mol_sys.run_md_sim(out_folder=os.path.join(out_folder, '03_prod'),
                           name='prod_' + sys_name,
                           mdp_template=equi_template_mdp,
                           mdp_options=free_ener_option_md, maxwarn=maxwarn,
                           monitor_tool=monitor_tool)
        pbar.update(prod_steps)

        return(os.path.join(out_folder, '03_prod/',
                            'prod_' + sys_name + '.xvg'))

    def compute_add_intermol_from_traj(self, ref_coor=None,
                                       rec_group='Protein', k=41.84):
        """ Compute intermolecular restraint from the last GmxSys trajectory.
        Get a coor object
        Get distance and angles
        """

        if ref_coor is None:
            ref_coor = self.ref_coor

        if not ref_coor.endswith('.pdb'):
            ref_sys = GmxSys(name='ref', coor_file=ref_coor)
            ref_sys.tpr = ref_coor
            ref_sys.convert_trj(traj=False, pbc='none')
            ref_coor_pdb = ref_sys.coor_file
        else:
            ref_coor_pdb = ref_coor
        coor = pdb_manip.Coor(ref_coor_pdb)

        # Center and align traj:
        self.gmxsys.convert_trj(
            select=rec_group + '\n System', center='yes')
        self.gmxsys.convert_trj(
            select=rec_group + '\n System', fit='rot+trans', pbc='none')

        # Get RMSF for ligand:
        lig_rmsf = self.gmxsys.get_rmsf([self.mol_name])
        # Sort df by RMSF
        lig_rmsf = lig_rmsf.sort_values(by=['RMSF'])
        # Get the stablest atom which is not Hydrogen
        for i, row in lig_rmsf.iterrows():
            atom = coor.atom_dict[int(row['Atom'])]
            if not atom['name'].startswith('H'):
                lig_atom_list = [int(row['Atom'])]
                break
        # Get connected atoms:
        lig_coor = coor.select_part_dict({'res_name': [self.mol_name]})
        atom_coor = pdb_manip.Coor()

        # Add 2 close atom consecutively
        for _ in range(2):
            atom_coor.atom_dict = {0: coor.atom_dict[lig_atom_list[-1]]}
            close_lig_atoms = lig_coor.get_index_dist_between(
                atom_coor, cutoff_min=1.0, cutoff_max=2.0)
            for i in close_lig_atoms:
                atom = coor.atom_dict[i]
                if not atom['name'].startswith('H') and i not in lig_atom_list:
                    lig_atom_list.append(i)
                    break

        # Atom index need to be +1 to be in gromacs numbering
        lig_atom_list = [index + 1 for index in lig_atom_list]
        logger.debug(f'Ligand atom indexes : {lig_atom_list}')

        # Get protein RMSF
        prot_rmsf = self.gmxsys.get_rmsf([rec_group])

        # Get backbone protein atom around the ligand:
        backbone = coor.select_part_dict({'name': ['N', 'C', 'O', 'CA']})
        # Create coor for the 3 ligand atoms
        atom_coor.atom_dict = {
            i: coor.atom_dict[lig_atom_list[i]] for i in range(3)}
        around_atom = backbone.get_index_dist_between(
            atom_coor, cutoff_min=1.0, cutoff_max=6.0)
        # Extract RMSF of contact atoms
        around_lig_df = prot_rmsf[prot_rmsf.Atom.isin(around_atom)]
        around_lig_df = around_lig_df.sort_values(by=['RMSF'])
        # Get residue of stablest atom:
        uniq_res = coor.atom_dict[
            around_lig_df.reset_index().loc[0, 'Atom']]['uniq_resid']
        rec_atom_list = coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': ['C']})
        rec_atom_list += coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': ['CA']})
        rec_atom_list += coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': ['N']})

        # Atom index need to be +1 to be in gromacs numbering
        rec_atom_list = [index + 1 for index in rec_atom_list]
        logger.debug(f'Receptor atom indexes : {rec_atom_list}')
        # Add C, CA, N

        self.add_intermol_restr_index(rec_atom_list, lig_atom_list,
                                      ref_coor_pdb, k=k)

    def add_intermol_restr_index(self, rec_index_list, lig_index_list,
                                 ref_coor, k=41.84, temp=300):
        """Compute the intermolecula restraints base on the
        self.coor file.

        :param rec_index_list: List of the receptor atom index
        :type rec_index_list: list

        :param lig_index_list: List of the ligand atom index
        :type lig_index_list: list

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.top_file

        Give three atoms for each receptor and ligand index list:
        R0, R1 R2 and L0 L1 L2
        Will define:
        - 1 bond:
            - R0-L0
        - 2 angles:
            - R0-L0-L1
            - R1-R0-L0
        - 2 dihedral angles:
            - R0-L0-L1-L2
            - R2-R1-R0-L0

        """
        bond_type = 6
        angle_type = 1
        dihed_type = 2

        # Get distance and angles

        if not ref_coor.endswith('.pdb'):
            ref_sys = GmxSys(name='ref', coor_file=ref_coor)
            ref_sys.tpr = ref_coor
            ref_sys.convert_trj(traj=False, pbc='none')
            coor = pdb_manip.Coor(ref_sys.coor_file)
        else:
            coor = pdb_manip.Coor(ref_coor)

        # R0-L0 (nm)
        # bond_df = self.get_dist([[rec_index_list[0], lig_index_list[0]]])
        # dist = bond_df.loc[1,:].mean() / 10 # Convert to nm

        # Index from gromacs starts at 1 and in pdb_manip at 0
        dist = pdb_manip.Coor.atom_dist(
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1]) / 10.0
        bond_list = [[rec_index_list[0],
                      lig_index_list[0],
                      bond_type,
                      round(dist, 3),
                      0,
                      round(dist, 3),
                      k * 100]]
        # angle_df = self.get_angle(
        #     [[rec_index_list[0], lig_index_list[0], lig_index_list[1]],
        #      [rec_index_list[1], rec_index_list[0], lig_index_list[0]]])
        # angle_1 = angle_df.loc[1,:].mean()
        # angle_2 = angle_df.loc[2,:].mean()

        # R0-L0-L1
        angle_1 = pdb_manip.Coor.atom_angle(
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1],
            coor.atom_dict[lig_index_list[1] - 1])
        # R1-R0-L0
        angle_2 = pdb_manip.Coor.atom_angle(
            coor.atom_dict[rec_index_list[1] - 1],
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1])
        angle_list = [[rec_index_list[0],
                       lig_index_list[0],
                       lig_index_list[1],
                       angle_type,
                       round(angle_1, 3), 0,
                       round(angle_1, 3), k]]
        angle_list += [[rec_index_list[1],
                        rec_index_list[0],
                        lig_index_list[0],
                        angle_type,
                        round(angle_2, 3), 0,
                        round(angle_2, 3), k]]

        # dihed_df = self.get_angle([[rec_index_list[0],
        #                             lig_index_list[0],
        #                             lig_index_list[1],
        #                             lig_index_list[2]],
        #                            [rec_index_list[2],
        #                             rec_index_list[1],
        #                             rec_index_list[0],
        #                             lig_index_list[0]],
        #                            [rec_index_list[1],
        #                             rec_index_list[0],
        #                             lig_index_list[0],
        #                             lig_index_list[1]]])
        # dihed_1 = dihed_df.loc[1,:].mean()
        # dihed_2 = dihed_df.loc[2,:].mean()
        # dihed_3 = dihed_df.loc[3,:].mean()
        # R0-L0-L1-L2
        dihed_1 = pdb_manip.Coor.atom_dihed_angle(
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1],
            coor.atom_dict[lig_index_list[1] - 1],
            coor.atom_dict[lig_index_list[2] - 1])
        # R2
        dihed_2 = pdb_manip.Coor.atom_dihed_angle(
            coor.atom_dict[rec_index_list[2] - 1],
            coor.atom_dict[rec_index_list[1] - 1],
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1])
        # R1-R0-L0-L1
        dihed_3 = pdb_manip.Coor.atom_dihed_angle(
            coor.atom_dict[rec_index_list[1] - 1],
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1],
            coor.atom_dict[lig_index_list[1] - 1])
        dihed_list = [[rec_index_list[0],
                       lig_index_list[0],
                       lig_index_list[1],
                       lig_index_list[2],
                       dihed_type,
                       round(dihed_1, 3), 0,
                       round(dihed_1, 3), k]]
        dihed_list += [[rec_index_list[2],
                        rec_index_list[1],
                        rec_index_list[0],
                        lig_index_list[0],
                        dihed_type,
                        round(dihed_2, 3), 0,
                        round(dihed_2, 3), k]]

        dihed_list += [[rec_index_list[1],
                        rec_index_list[0],
                        lig_index_list[0],
                        lig_index_list[1],
                        dihed_type,
                        round(dihed_3, 3), 0,
                        round(dihed_3, 3), k]]

        top = TopSys(self.gmxsys.top_file)

        top.add_intermolecular_restr(bond_list=bond_list,
                                     angle_list=angle_list,
                                     dihed_list=dihed_list)

        top.write_file(self.gmxsys.top_file[:-4] + '_rest.top')
        self.gmxsys.top_file = self.gmxsys.top_file[:-4] + '_rest.top'

        # Compute DG restr in water:
        V0 = 1660 * 1e-3  # nm³
        GK = KB * 1e-3 / 4.184  # Gas constant kcal/mol/K
        # K need to converted to Kcal (currently in KJ)
        k_bond = 1e2 * k / 4.184  # k*1e2 kcal/(nm*mol)
        k_angle = k / 4.184

        numerator = 8 * math.pi**2 * V0 * (k_bond * k_angle**5)**0.5
        denominator_1 = (dist**2 * math.sin(math.radians(angle_1))
                         * math.sin(math.radians(angle_2)))
        denominator_2 = (2 * math.pi * GK * temp)**3

        dg_rest_water = GK * temp * math.log(
            numerator / (denominator_1 * denominator_2))

        # Convert to kT
        dg_rest_water /= 0.593
        # Convert to self.unit
        dg_rest_water *= self.conv_fac

        logger.info(f'DG water restr = {dg_rest_water:.2f} {self.unit_name}')

        self.ener_rest_water = dg_rest_water

    def plot_intermol_restr(self, graph_out=None):
        """
        """
        fig, axs = plt.subplots(3, figsize=(8, 10))
        fig.suptitle(r'Inter Molecular Retraints')

        top = TopSys(self.gmxsys.top_file)

        dist_index_list = []
        for bond in top.inter_bond_list:
            dist = float(bond['rA'])
            dist_index_list.append([int(bond['ai']), int(bond['aj'])])
            axs[0].axvline(x=dist, c='0.5', linestyle='--')

        bond_df = self.gmxsys.get_dist(dist_index_list)

        min_dist = round(bond_df.iloc[:, 1].min(), 2)
        max_dist = round(bond_df.iloc[:, 1].max(), 2)
        print(min_dist, max_dist)
        # bins size of 0.1 Å or 0.01 nm
        bins_num = int((max_dist - min_dist) / 0.01)

        bond_df.iloc[:, 1:].plot.kde(ax=axs[0])
        bond_df.iloc[:, 1:].plot.hist(density=True, ax=axs[0], legend=False,
                                      bins=bins_num,
                                      range=(min_dist, max_dist))
        axs[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axs[0].set_xlabel('distance ($nm$)')
        axs[0].xaxis.set_label_coords(0.5, -0.1)

        angle_index_list = []
        for angle in top.inter_angl_list:
            angle_val = float(angle['thA'])
            angle_index_list.append([int(angle['ai']),
                                     int(angle['aj']),
                                     int(angle['ak'])])
            axs[1].axvline(x=angle_val, c='0.5', linestyle='--')

        angle_df = self.gmxsys.get_angle(angle_index_list)
        angle_df.iloc[:, 1:].plot.kde(ax=axs[1])
        angle_df.iloc[:, 1:].plot.hist(density=True, ax=axs[1],
                                       legend=False, alpha=0.5,
                                       bins=18 * 2, range=(0, 180))
        axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axs[1].set_xlabel(r'angle ($^{\circ}$)')
        axs[1].xaxis.set_label_coords(0.5, -0.1)

        dihe_index_list = []
        for angle in top.inter_dihe_list:
            dihe_val = float(angle['thA'])
            dihe_index_list.append([int(angle['ai']),
                                    int(angle['aj']),
                                    int(angle['ak']),
                                    int(angle['al'])])
            axs[2].axvline(x=dihe_val, c='0.5', linestyle='--')

        dihe_df = self.gmxsys.get_angle(dihe_index_list)
        dihe_df.iloc[:, 1:].plot.kde(ax=axs[2])
        dihe_df.iloc[:, 1:].plot.hist(density=True, ax=axs[2],
                                      legend=False, alpha=0.5,
                                      bins=36 * 2, range=(-180, 180))
        axs[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axs[2].set_xlabel(r'dihedral angle ($^{\circ}$)')
        axs[2].xaxis.set_label_coords(0.5, -0.1)

        if graph_out:
            plt.savefig(graph_out)

    @staticmethod
    def get_bar(xvg_file_list, bar_xvg='bar.xvg',
                barint_xvg='barint.xvg', hist_xvg='histogram.xvg',
                begin_time=0, end_time=-1,
                check_file_out=True, keep_ener_file=False):
        """Get energy of a system using ``gmx bar``.

        I don't know how to compute std like in gmx bar.
        """

        logger.info("- Extract bar energy")

        cmd_convert = os_command.Command([GMX_BIN, "bar",
                                          "-f", *xvg_file_list,
                                          "-o", bar_xvg,
                                          "-oi", barint_xvg,
                                          "-oh", hist_xvg,
                                          "-b", str(begin_time),
                                          "-e", str(end_time)])

        cmd_convert.display()
        cmd_convert.run()

        ener_pd = monitor.read_xvg(bar_xvg)

        if not keep_ener_file:
            os_command.delete_file(bar_xvg)
            os_command.delete_file(barint_xvg)
            os_command.delete_file(hist_xvg)

        return ener_pd

    def get_free_ener(self, begin_time=0, end_time=-1):
        """ Show free energy calculation output

        NEED TO FIX STD !!
        """
        self.table = FreeEner.get_bar(self.xvg_file_list,
                                      begin_time=begin_time,
                                      end_time=end_time)

        tot_ener = self.table['DG (kT)'].sum() * self.conv_fac
        tot_ener_std = sum(self.table['+/-']**2)**0.5 * self.conv_fac

        self.ener = tot_ener
        self.ener_std = tot_ener_std

        logger.info('\nDDG Tot   = {:5.2f} +/- {:.2f} {}\n'.format(
            tot_ener, tot_ener_std, self.unit_name))

        lambda_bond_num = len(self.lambda_restr)
        lambda_coul_num = len(self.lambda_coul)

        if self.lambda_restr:
            bond_contrib = self.table['DG (kT)'][
                :lambda_bond_num].sum() * self.conv_fac
            bond_contrib_std = sum(
                self.table['+/-'][:lambda_bond_num]**2)**0.5 * self.conv_fac

            logger.info('DG Restr = {:5.2f} +/- {:.2f} {}'.format(
                bond_contrib, bond_contrib_std, self.unit_name))

        if self.lambda_coul:
            coulomb_contrib = self.table['DG (kT)'][
                lambda_bond_num:lambda_bond_num + lambda_coul_num].sum()
            coulomb_contrib *= self.conv_fac
            coulomb_contrib_std = sum(self.table['+/-'][
                lambda_bond_num:lambda_bond_num + lambda_coul_num]**2)**0.5
            coulomb_contrib_std *= self.conv_fac
            logger.info('DG Coul  = {:5.2f} +/- {:.2f} {}'.format(
                coulomb_contrib, coulomb_contrib_std, self.unit_name))

        if self.lambda_vdw:
            vdw_contrib = self.table['DG (kT)'][
                lambda_bond_num + lambda_coul_num:].sum() * self.conv_fac
            vdw_contrib_std = sum(
                self.table['+/-'][lambda_bond_num + lambda_coul_num:]**2
                )**0.5 * self.conv_fac
            logger.info('DG LJ    = {:5.2f} +/- {:.2f} {}'.format(
                vdw_contrib, vdw_contrib_std, self.unit_name))

        return(tot_ener, tot_ener_std)

    def plot_convergence(self, graph_out=None, dt=10):
        """
        """
        try:
            self.compute_convergence_alchemlyb(dt=dt)
        except ImportError:
            self.compute_convergence_gbar(dt=dt)

        self.plot_convergence_graph(graph_out=graph_out)

    def compute_convergence_alchemlyb(self, dt=10):
        """
        """

        from alchemlyb.parsing.gmx import extract_dHdl
        from alchemlyb.estimators import TI

        files = self.xvg_file_list
        restr_num = len(self.lambda_restr)
        coul_num = len(self.lambda_coul)

        dHdl_tot = pd.concat([extract_dHdl(xvg, T=self.temp) for xvg in files])
        dHdl_tot_no_index = dHdl_tot.reset_index()
        ti_tot_list = []
        ti_tot_sd_list = []
        time_list = np.arange(dt, self.prod_time, dt)
        index_colnames = ['time', 'coul-lambda', 'vdw-lambda', 'bonded-lambda']
        self.convergence_data = {'time': time_list}

        if self.lambda_coul:
            dHdl_coul = pd.concat([extract_dHdl(xvg, T=self.temp)
                                   for xvg in files[
                                   restr_num:restr_num + coul_num]])
            dHdl_coul_no_index = dHdl_coul.reset_index()
            ti_coul_list = []
            ti_coul_sd_list = []

        if self.lambda_vdw:
            dHdl_vdw = pd.concat([extract_dHdl(xvg, T=self.temp)
                                  for xvg in files[restr_num + coul_num:]])
            dHdl_vdw_no_index = dHdl_vdw.reset_index()
            ti_vdw_list = []
            ti_vdw_sd_list = []

        if self.lambda_restr:
            dHdl_restr = pd.concat(
                [extract_dHdl(xvg, T=self.temp) for xvg in files[:restr_num]])
            dHdl_restr_no_index = dHdl_restr.reset_index()
            ti_restr_list = []
            ti_restr_sd_list = []

        for time in time_list:
            # Tot
            dHdl_local = dHdl_tot_no_index[
                (dHdl_tot_no_index.time < time) &
                (dHdl_tot_no_index.time > time - dt)]
            dHdl_local = dHdl_local.set_index(index_colnames)
            ti = TI().fit(dHdl_local)
            ti_tot_list.append(-self.conv_fac * ti.delta_f_.iloc[0, -1])
            ti_tot_sd_list.append(self.conv_fac * ti.d_delta_f_.iloc[0, -1])

            # Restr
            if self.lambda_restr:
                dHdl_local = dHdl_restr_no_index[
                    (dHdl_restr_no_index.time < time) &
                    (dHdl_restr_no_index.time > time - dt)]
                dHdl_local = dHdl_local.set_index(index_colnames)
                ti = TI().fit(dHdl_local)
                ti_restr_list.append(
                    -1 * self.conv_fac * ti.delta_f_.iloc[0, -1])
                ti_restr_sd_list.append(
                    self.conv_fac * ti.d_delta_f_.iloc[0, -1])

            # Coul
            if self.lambda_coul:
                dHdl_local = dHdl_coul_no_index[
                    (dHdl_coul_no_index.time < time) &
                    (dHdl_coul_no_index.time > time - dt)]
                dHdl_local = dHdl_local.set_index(index_colnames)
                ti = TI().fit(dHdl_local)
                ti_coul_list.append(-self.conv_fac * ti.delta_f_.iloc[0, -1])
                ti_coul_sd_list.append(
                    self.conv_fac * ti.d_delta_f_.iloc[0, -1])

            # VDW
            if self.lambda_vdw:
                dHdl_local = dHdl_vdw_no_index[
                    (dHdl_vdw_no_index.time < time) &
                    (dHdl_vdw_no_index.time > time - dt)]
                dHdl_local = dHdl_local.set_index(index_colnames)
                ti = TI().fit(dHdl_local)
                ti_vdw_list.append(-self.conv_fac * ti.delta_f_.iloc[0, -1])
                ti_vdw_sd_list.append(
                    self.conv_fac * ti.d_delta_f_.iloc[0, -1])

        self.convergence_data['tot'] = [ti_tot_list, ti_tot_sd_list]
        if self.lambda_restr:
            self.convergence_data['restr'] = [ti_restr_list, ti_restr_sd_list]
        if self.lambda_coul:
            self.convergence_data['coul'] = [ti_coul_list, ti_coul_sd_list]
        if self.lambda_vdw:
            self.convergence_data['vdw'] = [ti_vdw_list, ti_vdw_sd_list]

    def compute_convergence_gbar(self, dt=10):
        """
        """

        restr_num = len(self.lambda_restr)
        coul_num = len(self.lambda_coul)

        tot_list = []
        tot_sd_list = []
        coul_list = []
        coul_sd_list = []
        vdw_list = []
        vdw_sd_list = []
        restr_list = []
        restr_sd_list = []

        time_list = np.arange(dt, self.prod_time, dt)
        self.convergence_data = {'time': time_list}

        for time in time_list:
            local_table = FreeEner.get_bar(
                self.xvg_file_list, begin_time=time - dt, end_time=time)
            tot_contrib = local_table['DG (kT)'].sum() * self.conv_fac
            tot_list.append(-tot_contrib)
            tot_std_contrib = sum(self.table['+/-']**2)**0.5 * self.conv_fac
            tot_sd_list.append(tot_std_contrib)

            if self.lambda_restr:
                restr_contrib = local_table['DG (kT)'][:restr_num].sum()
                restr_contrib *= self.conv_fac
                restr_std_contrib = sum(local_table['+/-'][:restr_num]**2)**0.5
                restr_std_contrib *= self.conv_fac
                tot_list.append(-restr_contrib)
                tot_sd_list.append(restr_std_contrib)

            if self.lambda_coul:
                coul_contrib = local_table['DG (kT)'][
                    restr_num:restr_num + coul_num].sum()
                coul_contrib *= self.conv_fac
                coul_std_contrib = sum(local_table['+/-'][
                    restr_num:restr_num + coul_num]**2)**0.5
                coul_std_contrib *= self.conv_fac
                coul_list.append(-coul_contrib)
                coul_sd_list.append(coul_std_contrib)

            if self.lambda_vdw:
                vdw_contrib = local_table['DG (kT)'][
                    restr_num + coul_num:].sum()
                vdw_contrib *= self.conv_fac
                vdw_std_contrib = sum(
                    local_table['+/-'][restr_num + coul_num:]**2)**0.5
                vdw_std_contrib *= self.conv_fac
                vdw_list.append(-vdw_contrib)
                vdw_sd_list.append(vdw_std_contrib)

        self.convergence_data['tot'] = [tot_list, tot_sd_list]
        if self.lambda_restr:
            self.convergence_data['restr'] = [restr_list, restr_sd_list]
        if self.lambda_coul:
            self.convergence_data['coul'] = [coul_list, coul_sd_list]
        if self.lambda_vdw:
            self.convergence_data['vdw'] = [vdw_list, vdw_sd_list]

    def plot_convergence_graph(self, graph_out=None):
        """
        """

        fig, axs = plt.subplots(2, figsize=(8, 8), sharex=True)
        fig.suptitle(r'Convergence of $\Delta G$')
        fig.subplots_adjust(hspace=0)

        time_list = self.convergence_data['time']

        if self.lambda_restr:
            ti_restr_list = self.convergence_data['restr'][0]
            ti_restr_sd_list = self.convergence_data['restr'][1]
            axs[0].errorbar(
                time_list, ti_restr_list, ti_restr_sd_list,
                label=r'$\Delta G_{restr}$')
        if self.lambda_coul:
            ti_coul_list = self.convergence_data['coul'][0]
            ti_coul_sd_list = self.convergence_data['coul'][1]
            axs[0].errorbar(
                time_list, ti_coul_list, ti_coul_sd_list,
                label=r'$\Delta G_{coul}$')
        if self.lambda_vdw:
            ti_vdw_list = self.convergence_data['vdw'][0]
            ti_vdw_sd_list = self.convergence_data['vdw'][1]
            axs[0].errorbar(
                time_list, ti_vdw_list, ti_vdw_sd_list,
                label=r'$\Delta G_{vdw}$')
        axs[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        ti_tot_list = self.convergence_data['tot'][0]
        ti_tot_sd_list = self.convergence_data['tot'][1]
        axs[1].errorbar(
            time_list, ti_tot_list, ti_tot_sd_list, label=r'$\Delta G_{tot}$')
        axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        axs[1].set_ylabel(r'$\Delta G \; ({})$'.format(self.unit_graph))
        axs[1].yaxis.set_label_coords(-0.1, 1)
        axs[1].set_xlabel(r'$time \; (ps)$')

        if graph_out:
            plt.savefig(graph_out)

    @staticmethod
    def symmetry_correction(smile, temp=300):
        r""" Compute symmetry correction
        $\Delta_{sym} = −k T ln(\sigma) $

        return value in kcal/mol


        >>> FreeEner.symmetry_correction('c1ccccc1') #doctest: +ELLIPSIS
        -1.4814...
        >>> FreeEner.symmetry_correction('c1ccccc1O') #doctest: +ELLIPSIS
        -0.4132...
        >>> FreeEner.symmetry_correction('c1cc(O)ccc1O') #doctest: +ELLIPSIS
        -0.8264...
        """

        try:
            from rdkit import Chem
            # from rdkit.Chem import rdmolfiles
        except ImportError:
            logger.warning('WARNING !!!! \n'
                           'Could not load rdkit \nInstall it using conda:\n'
                           'conda install -c conda-forge rdkit\n'
                           'Symmetry correction set to 0\n')
            return(0.0)

        x = Chem.MolFromSmiles(smile)
        # z = list(rdmolfiles.CanonicalRankAtoms(x, breakTies=False))
        matches = x.GetSubstructMatches(x, uniquify=False)
        sigma = len(matches)
        logger.warning('WARNING !!!! \n'
                       'symmetry_correction() is quite experimental\n'
                       'seems to work for small molecules, use it with'
                       ' caution\n')
        logger.info(f'For molecule {smile} Symmetry number'
                    ' sigma is set to {sigma}')

        GK = - KB * 1e-3 / 4.184  # Gas constant kcal/mol/K

        return GK * temp * math.log(sigma)
