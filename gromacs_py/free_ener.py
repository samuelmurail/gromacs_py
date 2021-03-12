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
from .gmx import GROMACS_MOD_DIRNAME
from .tools import monitor, ambertools

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Development"


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

    gmx.show_log()


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
    """ Free Energy encapsulation class.

    This class can be used to launch and analyze free energy
    calculations using Free Energy Perturbation.

    :param out_folder: output folder
    :type out_folder: str

    :param lambda_coul: lambda points for Coulomb
    :type lambda_coul: list of float

    :param lambda_vdw: lambda points for Lennard Jones
    :type lambda_vdw: list of float

    :param lambda_restr: lambda points for restraints
    :type lambda_restr: list of float

    :param xvg_file_list: List of free energy xvg files
    :type xvg_file_list: list of string

    :param lambda_sys_list: List of lambda GmxSys
    :type lambda_sys_list: list of string

    :param temp: Temperature
    :type temp: float

    :param smile: Ligand SMILE
    :type smile: str

    """

    def __init__(self, mol_name, out_folder, unit='kcal'):
        self.mol_name = mol_name[:3]
        if len(mol_name) > 4:
            logger.warning(
                f'Molecule name {mol_name} is longer than 3 characters,'
                f'name changed to {self.mol_name}')
        self.unit = unit
        self.out_folder = out_folder
        # Lambda
        self.lambda_coul = []
        self.lambda_vdw = []
        self.lambda_restr = []
        # Xvg file list
        self.xvg_file_list = []
        self.lambda_sys_list = []
        self.temp = None
        self.smile = None

    @property
    def conv_fac(self):
        """ Conversion factor from kT to
        self.unit
        """
        return(FreeEner.get_conv_fac(self.unit, self.temp))

    @staticmethod
    def get_conv_fac(unit, temp):
        """ Conversion factor from kT to
        self.unit
        """
        if unit == 'kcal':
            return(0.593)
        elif unit == 'kJ':
            return(4.184 * 0.593)
        elif unit == 'kT':
            return(1)
        elif unit == 'logP':
            return(-1.365679 * KB * temp / 1000)

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
        """ Return unit as math latex for matplotlib
        label purpose.

        """
        if self.unit == 'kcal':
            return(r'kcal \cdot mol^{-1}')
        elif self.unit == 'kJ':
            return(r'kJ \cdot mol^{-1}')
        elif self.unit == 'kT':
            return(r'k_{B}T')
        elif self.unit == 'logP':
            return(r'log \; P')

    def water_box_from_SMILE(self, smile, method_3d='rdkit',
                             iter_num=5000, box_dist=1.1):
        """ Create a water box coordinates and toplogie with a
        molecule defined by its SMILE.

        :param smile: Molecule's SMILE
        :type smile: str

        :param method_3d: Method to compute 3D coordinates
        :type method_3d: str, default='rdkit'

        :param iter_num: Iteration step number for 3D coordinate
            computation
        :type iter_num: int, default=5000

        :param box_dist: Ditance to egde box (nm)
        :type box_dist: float, default=1.1 nm

        .. Note::
            Default box dist 1.1 nm was taken as minimal distance for
            CH4 molecule, to allow domain decomposition with `gmx mdrun`.

        """

        self.smile = smile

        os.makedirs(self.out_folder, exist_ok=True)

        # Create 3d structure from SMILE
        pdb_file = os.path.join(self.out_folder, f'{self.mol_name}.pdb')
        charge = ambertools.smile_to_pdb(smile, pdb_file,
                                         method_3d=method_3d,
                                         mol_name=self.mol_name,
                                         iter_num=iter_num)
        logger.info(f'ligand charge is {charge}')

        # Topologie and system creation
        self.gmxsys = gmx.GmxSys(name=self.mol_name, coor_file=pdb_file)
        self.gmxsys.prepare_top_ligand(
            out_folder=os.path.join(self.out_folder, 'mol_top'),
            ff='amber99sb-ildn', include_mol={self.mol_name: smile})

        self.gmxsys.create_box(name=None, dist=box_dist)
        self.gmxsys.solvate_box(
            out_folder=os.path.join(self.out_folder, 'mol_water_top'))

    def octanol_box_from_SMILE(self, smile, method_3d='rdkit',
                               iter_num=5000, box_dist=1.3):
        """ Create an octanol box coordinates and toplogie with a
        molecule defined by its SMILE.

        :param smile: Molecule's SMILE
        :type smile: str

        :param method_3d: Method to compute 3D coordinates
        :type method_3d: str, default='rdkit'

        :param iter_num: Iteration step number for 3D coordinate
            computation
        :type iter_num: int, default=5000

        :param box_dist: Ditance to egde box (nm)
        :type box_dist: float, default=1.3 nm

        .. Note::
            Default box dist 1.3 nm was taken as minimal distance for
            CH4 molecule, to allow domain decomposition with `gmx mdrun`.
        """

        self.smile = smile

        os.makedirs(self.out_folder, exist_ok=True)

        # Create 3d structure from SMILE
        pdb_file = os.path.join(self.out_folder, f'{self.mol_name}.pdb')
        charge = ambertools.smile_to_pdb(smile, pdb_file,
                                         method_3d=method_3d,
                                         mol_name=self.mol_name,
                                         iter_num=iter_num)
        logger.info(f'ligand charge is {charge}')

        # Topologie and system creation
        self.gmxsys = gmx.GmxSys(name=self.mol_name, coor_file=pdb_file)
        self.gmxsys.prepare_top_ligand(
            out_folder=os.path.join(self.out_folder, 'mol_top'),
            ff='amber99sb-ildn', include_mol={self.mol_name: smile})

        self.gmxsys.create_box(name=None, dist=box_dist)
        self.gmxsys.solvate_box(
            out_folder=os.path.join(self.out_folder, 'mol_octanol_top'),
            cs=os.path.join(GROMACS_MOD_DIRNAME,
                            "template/octanol/prod_OCO.gro"))

        # Add octanol topologie:
        sys_top = gmx.TopSys(self.gmxsys.top_file)
        atomtypes_itp = os.path.join(GROMACS_MOD_DIRNAME,
                                     "template/octanol/OCO_GMX_atomtypes.itp")
        sys_top.add_atomtypes(atomtypes_itp)

        mol_itp = os.path.join(GROMACS_MOD_DIRNAME,
                               "template/octanol/OCO_GMX.itp")
        sys_top.add_mol_itp(mol_itp)

        sys_top.write_file(self.gmxsys.top_file)

    def equilibrate_solvent_box(self, em_steps=10000, dt=0.002, prod_time=10.0,
                                short_steps=50000, temp=300.0):
        """ Equilibrate a solvent (water, octanol) box.

        :param em_steps: Energy minimisation step number
        :type em_steps: int, default=10000

        :param dt: Integration time step (ps)
        :type dt: float, default=0.002 ps

        :param prod_time: Equilibration time (ns)
        :type prod_time: float, default=10.0 ns

        :param short_steps: Short equilibration steps
        :type short_steps: int, default=50000

        :param temp: Temperature of equilibration
        :type temp: float, default=300.0

        """

        # EM
        self.gmxsys.em_2_steps(out_folder=os.path.join(
                                self.out_folder, 'mol_solv_em'),
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
            out_folder=os.path.join(self.out_folder, 'mol_solv_prod'),
            name='short_equi',
            nsteps=short_steps, dt=0.0005,
            maxwarn=1, **mdp_options)
        # second step
        prod_steps = 1000 * prod_time / dt
        self.gmxsys.production(
            out_folder=os.path.join(self.out_folder, 'mol_solv_prod'),
            nsteps=prod_steps,
            dt=dt, maxwarn=1, **mdp_options)

    def prepare_complex_pdb(self, pdb_in, smile, ff='amber99sb-ildn'):
        """ Prepare topologie from a pdb file, create box around and solvate it.

        :param pdb_in: Input PDB file
        :type pdb_in: str

        :param smile: Ligand's SMILE
        :type smile: str

        :param ff: Forcefield
        :type ff: str, default='amber99sb-ildn'

        """

        self.gmxsys = gmx.GmxSys(name=self.mol_name, coor_file=pdb_in)
        self.gmxsys.prepare_top(out_folder=os.path.join(
                                    self.out_folder, 'complex_top'),
                                ff=ff, include_mol={self.mol_name: smile})
        self.gmxsys.solvate_add_ions(out_folder=os.path.join(
                                        self.out_folder, 'sys_top'),
                                     ion_C=0.15, maxwarn=3,
                                     create_box_flag=True)

    def equilibrate_complex(self, em_steps=10000, HA_time=0.25,
                            CA_time=0.5, CA_LOW_time=1.0,
                            dt=0.002, dt_HA=0.001, temp=300,
                            receptor_grp='Protein',
                            short_steps=50000):
        """ Equilibrate a receptor-ligand complex system.

        :param em_steps: Energy minimisation step number, default=10000
        :type em_steps: int

        :param HA_time: Heavy atoms restraints equilibration time (ns)
        :type HA_time: float, default=0.25

        :param CA_time: Alpha carbon atoms restraints equilibration
            time (ns)
        :type CA_time: float, default=0.5

        :param CA_LOW_time: Low alpha carbon atoms restraints equilibration
            time (ns)
        :type CA_LOW_time: float, default=1.0

        :param dt: Integration time step (ps)
        :type dt: float, default=0.002 ps

        :param dt_HA: Integration time step (ps)
        :type dt_HA: float, default=0.001 ps

        :param temp: Temperature of equilibration (K)
        :type temp: float, default=300.0 K

        :param receptor_grp: Receptor group (for temperature coupling)
        :type receptor_grp: str, default='Protein'

        :param short_steps: Short equilibration steps
        :type short_steps: int, default=50000

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
        mdp_options = {'nsteps': short_steps,
                       'define': '-DPOSRES', 'dt': 0.0005,
                       'tc-grps': '{} {} Water_and_ions'.format(
                            receptor_grp,
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

        mdp_options = {}
        if receptor_grp != 'Protein':
            index_dict = self.gmxsys.get_index_dict()
            self.gmxsys.add_ndx('{}|{} \nq\n'.format(
                    index_dict[receptor_grp],
                    index_dict[self.mol_name]),
                check_file_out=False)

            mdp_options = {'tc-grps': '{}_{} Water_and_ions'.format(
                            receptor_grp,
                            self.mol_name)}
        print(mdp_options)

        self.gmxsys.equi_three_step(out_folder=os.path.join(
                                     self.out_folder, 'sys_equi'),
                                    nsteps_HA=1000 * HA_time / dt_HA,
                                    nsteps_CA=1000 * CA_time / dt,
                                    nsteps_CA_LOW=1000 * CA_LOW_time / dt,
                                    dt=dt, dt_HA=dt_HA, maxwarn=3,
                                    pdb_restr=self.ref_coor,
                                    **mdp_options)

    def run(self, lambda_coul_list, lambda_vdw_list, lambda_restr_list=[],
            mbar=False, dir_name='free_ener_run',
            em_steps=5000, nvt_time=10, npt_time=10, prod_time=100,
            dt=0.002, temp=300.0,
            temp_groups='Protein non-Protein', maxwarn=1,
            monitor_tool=monitor.PROGRESS_BAR):
        """Compute free energy to uncouple a molecule to a
        system.

        :param dir_name: path of the output folder
        :type dir_name: str, default='free_ener_run'

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

            * self.gmxsys.coor_file
            * self.gmxsys.top_file
            * self.gmxsys.nt
            * self.gmxsys.ntmpi
            * self.gmxsys.gpu_id

        **Object field(s) changed:**

            * self.lambda_sys_list
            * self.lambda_coul
            * self.lambda_vdw
            * self.lambda_restr
            * self.prod_time
            * self.xvg_file_list

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
                               'lincs_order': 8,  # Try to fix the LINCS Errors
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
                i + 1, lambda_restr_num + lambda_coul_num + lambda_vdw_num))

            lambda_sys = FreeEner.compute_lambda_point(
                self.gmxsys, i,
                self.mol_name,
                os.path.join(self.out_folder,
                             dir_name),
                free_ener_option_md, pbar,
                mbar=mbar,
                em_steps=em_steps,
                nvt_steps=nvt_steps,
                npt_steps=npt_steps,
                prod_steps=prod_steps,
                maxwarn=maxwarn,
                monitor_tool=monitor_tool)
            self.lambda_sys_list.append(lambda_sys)
            self.xvg_file_list.append(lambda_sys.xvg)

        self.temp = temp
        self.lambda_coul = list(lambda_coul_list)
        self.lambda_vdw = list(lambda_vdw_list)
        self.lambda_restr = list(lambda_restr_list)
        self.prod_time = prod_time

    @staticmethod
    def compute_lambda_point(gmx_sys, lambda_iter, mol_name, out_folder,
                             free_ener_option, pbar, mbar,
                             em_steps, nvt_steps, npt_steps,
                             prod_steps, maxwarn=1,
                             monitor_tool=monitor.PROGRESS_BAR):
        """ Run the different steps of a single lambda point.

        :param gmx_sys: Gmx System to start from
        :type gmx_sys: GmxSys

        :param lambda_iter: lambda point
        :type lambda_iter: int

        :param mol_name: Molecule residue name
        :type mol_name: str

        :param out_folder: path of the output folder
        :type out_folder: str

        :param free_ener_option: Mdp options
        :type free_ener_option: dict

        :param pbar: progress bar object
        :type pbar: tqmd bar

        :param mbar: MBAR flag
        :type mbar: bool

        :param em_steps: number of minimisation steps
        :type em_steps: int

        :param nvt_steps: number of NVT equilibration steps
        :type nvt_steps: int

        :param npt_steps:  number of NPT equilibration steps
        :type npt_steps: int

        :param prod_steps: number of production steps
        :type prod_steps: int

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type rerun: dict, default=None

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

        return(mol_sys)

    def extend_lambda_prod(self, prod_time):
        """ Extend free energy production.

        :param prod_time: production time (ns)
        :type prod_time: float

        """

        dt = float(self.lambda_sys_list[0].get_mdp_dict()['dt'])
        nsteps = prod_time / dt
        lambda_num = len(self.lambda_restr) + len(self.lambda_coul) +\
            len(self.lambda_vdw)
        xvg_file_list = []

        for i in range(lambda_num):
            logger.info(f'Compute lambda {i+1} / {lambda_num}')

            self.lambda_sys_list[i].extend_sim(nsteps=nsteps)
            output = self.lambda_sys_list[i].get_all_output()
            xvg_file_list += output['xvg']

        self.prod_time = prod_time
        self.xvg_file_list = xvg_file_list

    def align_ref_traj(self, rec_group='Protein'):
        """
        """

        # Center and align traj:
        self.gmxsys.convert_trj(
            select=rec_group + '\n System', center='yes')
        self.gmxsys.convert_trj(
            select=rec_group + '\n System', fit='rot+trans', pbc='none')

    def compute_add_intermol_from_traj(self, ref_coor=None,
                                       rec_group='Protein', k=41.84,
                                       cutoff_prot=6.0):
        """ Compute intermolecular restraint from the last GmxSys trajectory.
        Get a coor object
        Get distance and angles
        """

        if ref_coor is None:
            ref_coor = self.ref_coor

        self.align_ref_traj(rec_group=rec_group)

        lig_atom_list = self.get_ligand_atoms(ref_coor)

        rec_atom_list = self.get_protein_atoms_from_rmsf(
            ref_coor,
            lig_atom_list,
            rec_group=rec_group,
            cutoff_max=cutoff_prot)

        self.add_intermol_restr_index(rec_atom_list, lig_atom_list,
                                      ref_coor, k=k)

    def get_protein_atoms_from_rmsf(self, ref_coor, lig_atom_list,
                                    rec_group='Protein',
                                    cutoff_max=6.0):
        """
        """
        coor = pdb_manip.Coor(ref_coor)

        # Get protein RMSF
        prot_rmsf = self.gmxsys.get_rmsf([rec_group])

        # Get backbone protein atom around the ligand:
        if rec_group == 'Protein':
            bk_atoms_list = ['C', 'CA', 'N', 'O']
        elif rec_group == 'DNA':
            bk_atoms_list = ['C4\'', 'C3\'', 'O3\'', 'O5\'', 'P', 'C5\'']
        backbone = coor.select_part_dict({'name': bk_atoms_list})
        # Create coor for the 3 ligand atoms
        atom_coor = pdb_manip.Coor()
        atom_coor.atom_dict = {
            i: coor.atom_dict[lig_atom_list[i] - 1] for i in range(3)}
        around_atom = backbone.get_index_dist_between(
            atom_coor, cutoff_min=0.0, cutoff_max=cutoff_max)
        # Extract RMSF of contact atoms
        around_lig_df = prot_rmsf[prot_rmsf.Atom.isin(around_atom)]
        around_lig_df = around_lig_df.sort_values(by=['RMSF'])
        # Get residue of stablest atom:
        uniq_res = coor.atom_dict[
            around_lig_df.reset_index().loc[0, 'Atom']]['uniq_resid']
        rec_atom_list = coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': [bk_atoms_list[0]]})
        rec_atom_list += coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': [bk_atoms_list[1]]})
        rec_atom_list += coor.get_index_selection(
            {'uniq_resid': [uniq_res], 'name': [bk_atoms_list[2]]})

        # Atom index need to be +1 to be in gromacs numbering
        rec_atom_list = [index + 1 for index in rec_atom_list]
        logger.debug(f'Receptor atom indexes : {rec_atom_list}')

        return(rec_atom_list)

    def get_protein_atoms_from_res(self, resid,
                                   rec_group='Protein'):
        """
        """

        coor = pdb_manip.Coor(self.ref_coor)

        # Get backbone protein atom around the ligand:
        if rec_group == 'Protein':
            bk_atoms_list = ['C', 'CA', 'N']
        elif rec_group == 'DNA':
            bk_atoms_list = ['C4\'', 'C3\'', 'O3\'']

        rec_atom_list = coor.get_index_selection(
            {'res_num': [resid], 'name': [bk_atoms_list[0]]})
        rec_atom_list += coor.get_index_selection(
            {'res_num': [resid], 'name': [bk_atoms_list[1]]})
        rec_atom_list += coor.get_index_selection(
            {'res_num': [resid], 'name': [bk_atoms_list[2]]})

        # Atom index need to be +1 to be in gromacs numbering
        rec_atom_list = [index + 1 for index in rec_atom_list]
        logger.debug(f'Receptor atom indexes : {rec_atom_list}')

        return(rec_atom_list)

    def get_ligand_atoms(self, ref_coor):
        """
        """

        coor = pdb_manip.Coor(ref_coor)
        # Get RMSF for ligand:
        lig_rmsf = self.gmxsys.get_rmsf([self.mol_name])
        # Sort df by RMSF
        lig_rmsf = lig_rmsf.sort_values(by=['RMSF'])
        # Get the stablest atom which is not Hydrogen
        for i, row in lig_rmsf.iterrows():
            atom = coor.atom_dict[int(row['Atom']) - 1]
            if not atom['name'].startswith('H'):
                lig_atom_list = [int(row['Atom']) - 1]
                break
        # Get connected atoms:
        lig_coor = coor.select_part_dict({'res_name': [self.mol_name]})
        atom_coor = pdb_manip.Coor()

        # Add 2 close atom consecutively
        cutoff_max = 1.5
        while len(lig_atom_list) < 3:
            cutoff_max += 0.1
            atom_coor.atom_dict = {0: coor.atom_dict[lig_atom_list[-1]]}
            close_lig_atoms = lig_coor.get_index_dist_between(
                atom_coor, cutoff_min=0.1, cutoff_max=cutoff_max)
            if len(close_lig_atoms) > 0:
                for _, row in lig_rmsf.iterrows():
                    num = int(row['Atom']) - 1
                    atom = coor.atom_dict[num]
                    if num in close_lig_atoms and not (
                            atom['name'].startswith('H')
                            or num in lig_atom_list):
                        lig_atom_list.append(num)
                        break

        # Atom index need to be +1 to be in gromacs numbering
        lig_atom_list = [index + 1 for index in lig_atom_list]
        logger.debug(f'Ligand atom indexes : {lig_atom_list}')
        print(f'Ligand atom indexes : {lig_atom_list}')

        # In rare case there is more than 3 atoms
        return(lig_atom_list[:3])

    def add_intermol_restr_index(self, rec_index_list, lig_index_list,
                                 ref_coor, k=41.84, temp=300):
        r"""Compute and add the intermolecular restraints based on the
        ref_coor distance and angles.

        Give three atoms for each receptor and ligand index list:
        :math:`R_0`, :math:`R_1`, :math:`R_2` and :math:`L_0`, :math:`L_1`,
        :math:`L_2`
        Will define:

        * 1 bond:

            * :math:`R_0 - L_0`

        * 2 angles:

            * :math:`R_0 - L_0 - L_1`
            * :math:`R_1 - R_0 - L_0`

        * 2 dihedral angles:

            * :math:`R_0 - L_0 - L_1 - L_2`
            * :math:`R_2 - R_1 - R_0 - L_0`

        :param rec_index_list: List of the receptor atom index
        :type rec_index_list: list

        :param lig_index_list: List of the ligand atom index
        :type lig_index_list: list

        :param ref_coor: Reference coordinates file
        :type ref_coor: str

        :param k: intermolecular force constant,
            (default= :math:`41.84 \; kcal \, mol^{-1} \, nm^{-2}`)
        :type k: float

        :param temp: Temperature defalult=300 K
        :type temp: float

        **Object requirement(s):**

            * self.gmxsys.coor_file
            * self.gmxsys.top_file

        **Object field(s) changed:**

            * self.gmxsys.top_file

        """
        bond_type = 6
        angle_type = 1
        dihed_type = 2

        coor = pdb_manip.Coor(ref_coor)

        # R0-L0
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

        # R0-L0-L1-L2
        dihed_1 = pdb_manip.Coor.atom_dihed_angle(
            coor.atom_dict[rec_index_list[0] - 1],
            coor.atom_dict[lig_index_list[0] - 1],
            coor.atom_dict[lig_index_list[1] - 1],
            coor.atom_dict[lig_index_list[2] - 1])
        # R2-R1-R0-L0
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

        top = gmx.TopSys(self.gmxsys.top_file)

        top.add_intermolecular_restr(bond_list=bond_list,
                                     angle_list=angle_list,
                                     dihed_list=dihed_list)

        top.write_file(self.gmxsys.top_file[:-4] + '_rest.top')
        self.gmxsys.top_file = self.gmxsys.top_file[:-4] + '_rest.top'

    def get_water_restr(self, temp=300):
        r""" Compute ligand restaint energy in water using
        Boresh et al. equation:

        :math:`\Delta G_{restr} = RT \ln \left( \dfrac{8 \pi^2 V^0}
        {r^2_0 \sin{\theta_{a 0}} \sin{\theta_{b 0}}} \dfrac{
        \sqrt{k_r k_{\theta_a} k_{\theta_b} k_{\tau_\alpha} k_{\tau_\beta}
        k_{\tau_\gamma}}} {2 \pi KT^3} \right)`

        """

        top = gmx.TopSys(self.gmxsys.top_file)

        dist = float(top.inter_bond_list[0]['rA'])
        k_bond = float(top.inter_bond_list[0]['kB'])

        k_angle = float(top.inter_angl_list[0]['kB'])
        angle_1 = float(top.inter_angl_list[0]['thA'])
        angle_2 = float(top.inter_angl_list[1]['thA'])

        # Compute DG restr in water:
        V0 = 1660 * 1e-3  # nm³
        GK = KB * 1e-3 / 4.184  # Gas constant kcal/mol/K
        # K need to converted to Kcal (currently in KJ)
        k_bond /= 4.184  # k*1e2 kcal/(nm*mol)
        k_angle /= 4.184

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

        return(dg_rest_water)

    def show_intermol_restr(self):
        """ Show traj with atom implied in intermolecular
        restraints using nglview library.
        """

        view = self.gmxsys.view_traj()
        view.clear_representations()
        view.add_representation("cartoon")

        # Get distance index
        top = gmx.TopSys(self.gmxsys.top_file)

        atom_list = []
        for bond in top.inter_bond_list:
            atom_a = int(bond['ai']) - 1
            atom_b = int(bond['aj']) - 1
            view.add_representation("distance",
                                    atomPair=[(f'@{atom_a}',
                                               f'@{atom_b}')])
            atom_list += [atom_a, atom_b]

        # Get all atoms:
        for angle in top.inter_dihe_list:
            atom_a = int(angle['ai']) - 1
            atom_b = int(angle['aj']) - 1
            atom_c = int(angle['ak']) - 1
            atom_d = int(angle['al']) - 1
            atom_list += [atom_a, atom_b, atom_c, atom_d]
        # Get unique value:
        atom_list = list(set(atom_list))

        # Get residue number:
        coor = pdb_manip.Coor(self.ref_coor)
        residue_list = []
        for atom in atom_list:
            residue_list.append(coor.atom_dict[atom]['res_num'])
        residue_list = list(set(residue_list))

        selection = f'[{self.mol_name}]'
        if len(residue_list) > 0:
            selection += ' or '
            selection += ','.join([str(i) for i in residue_list])
        print('residue list:', residue_list)

        view.add_representation("licorice",
                                selection=selection)
        atom_list_sel = '@' + ','.join([str(i) for i in atom_list])
        view.add_representation("ball+stick",
                                selection=atom_list_sel,
                                aspectRatio=4.0)

        return(view)

    def plot_intermol_restr(self, graph_out=None):
        """
        """
        fig, axs = plt.subplots(3, figsize=(8, 10))
        fig.suptitle(r'Inter Molecular Retraints')

        width_ratio = 0.7

        top = gmx.TopSys(self.gmxsys.top_file)

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
        axs[0].grid()
        # To avoid legend being cut:
        box = axs[0].get_position()
        axs[0].set_position([box.x0, box.y0,
                             box.width * width_ratio,
                             box.height])

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
        axs[1].set_xlim(0, 180)
        axs[1].grid()
        # To avoid legend being cut:
        box = axs[1].get_position()
        axs[1].set_position([box.x0, box.y0,
                             box.width * width_ratio,
                             box.height])

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
        axs[2].set_xlim(-180, 180)
        axs[2].grid()
        # To avoid legend being cut:
        box = axs[2].get_position()
        axs[2].set_position([box.x0, box.y0,
                             box.width * width_ratio,
                             box.height])

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

    def get_free_ener(self, begin_time=0, end_time=-1, unit=None):
        """ Show free energy calculation output

        NEED TO FIX STD !!
        """
        if unit is None:
            conv_fac = self.conv_fac
        else:
            conv_fac = FreeEner.get_conv_fac(unit, self.temp)

        self.table = FreeEner.get_bar(self.xvg_file_list,
                                      begin_time=begin_time,
                                      end_time=end_time)

        tot_ener = self.table['DG (kT)'].sum() * conv_fac
        tot_ener_std = sum(self.table['+/-']**2)**0.5 * conv_fac

        self.ener = tot_ener
        self.ener_std = tot_ener_std

        logger.info('\nDDG Tot   = {:5.2f} +/- {:.2f} {}\n'.format(
            tot_ener, tot_ener_std, self.unit_name))

        lambda_bond_num = len(self.lambda_restr)
        lambda_coul_num = len(self.lambda_coul)

        if self.lambda_restr:
            bond_contrib = self.table['DG (kT)'][
                :lambda_bond_num].sum() * conv_fac
            bond_contrib_std = sum(
                self.table['+/-'][:lambda_bond_num]**2)**0.5 * conv_fac

            logger.info('DG Restr = {:5.2f} +/- {:.2f} {}'.format(
                bond_contrib, bond_contrib_std, self.unit_name))

        if self.lambda_coul:
            coulomb_contrib = self.table['DG (kT)'][
                lambda_bond_num:lambda_bond_num + lambda_coul_num].sum()
            coulomb_contrib *= conv_fac
            coulomb_contrib_std = sum(self.table['+/-'][
                lambda_bond_num:lambda_bond_num + lambda_coul_num]**2)**0.5
            coulomb_contrib_std *= conv_fac
            logger.info('DG Coul  = {:5.2f} +/- {:.2f} {}'.format(
                coulomb_contrib, coulomb_contrib_std, self.unit_name))

        if self.lambda_vdw:
            vdw_contrib = self.table['DG (kT)'][
                lambda_bond_num + lambda_coul_num:].sum() * conv_fac
            vdw_contrib_std = sum(
                self.table['+/-'][lambda_bond_num + lambda_coul_num:]**2
                )**0.5 * conv_fac
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

        return(self.plot_convergence_graph(graph_out=graph_out))

    def compute_convergence_alchemlyb(self, dt=10):
        """
        """

        from alchemlyb.parsing.gmx import extract_dHdl
        from alchemlyb.estimators import TI

        files = self.xvg_file_list
        # NEED TO SORTS FILES IN CASE OF RESTARTS
        restr_num = len(self.lambda_restr)
        coul_num = len(self.lambda_coul)
        vdw_num = len(self.lambda_vdw)

        dHdl_tot = pd.concat([extract_dHdl(xvg, T=self.temp) for xvg in files])
        dHdl_tot_no_index = dHdl_tot.reset_index()
        # Remove duplicates due to restarts
        dHdl_tot_no_index = dHdl_tot_no_index.drop_duplicates(
            subset=['bonded-lambda', 'coul-lambda', 'vdw-lambda', 'time'],
            keep='last')
        # Sort values
        dHdl_tot_no_index = dHdl_tot_no_index.sort_values(
            by=['bonded-lambda', 'coul-lambda', 'vdw-lambda', 'time'])
        # Reset index again:
        dHdl_tot_no_index = dHdl_tot_no_index.reset_index(drop=True)
        # return(dHdl_tot_no_index)

        ti_tot_list = []
        ti_tot_sd_list = []
        time_list = np.arange(dt, self.prod_time + dt, dt)
        index_colnames = ['time', 'coul-lambda', 'vdw-lambda', 'bonded-lambda']
        self.convergence_data = {'time': time_list}

        sim_num_val = len(dHdl_tot_no_index) / (restr_num + coul_num + vdw_num)

        if self.lambda_coul:
            dHdl_coul_no_index = dHdl_tot_no_index[
                max(0, int(sim_num_val * (restr_num - 1))):
                int(sim_num_val * (restr_num + coul_num))]
            ti_coul_list = []
            ti_coul_sd_list = []

        if self.lambda_vdw:
            dHdl_vdw_no_index = dHdl_tot_no_index[
                int(sim_num_val * (restr_num + coul_num - 1)):]
            ti_vdw_list = []
            ti_vdw_sd_list = []

        if self.lambda_restr:
            dHdl_restr_no_index = dHdl_tot_no_index[
                :int(sim_num_val * restr_num)]
            ti_restr_list = []
            ti_restr_sd_list = []

        for time in time_list:
            # Tot
            dHdl_local = dHdl_tot_no_index[
                (dHdl_tot_no_index.time < time) &
                (dHdl_tot_no_index.time >= time - dt)]
            dHdl_local = dHdl_local.set_index(index_colnames)
            ti = TI().fit(dHdl_local)
            ti_tot_list.append(-self.conv_fac * ti.delta_f_.iloc[0, -1])
            ti_tot_sd_list.append(self.conv_fac * ti.d_delta_f_.iloc[0, -1])

            # Restr
            if self.lambda_restr:
                dHdl_local = dHdl_restr_no_index[
                    (dHdl_restr_no_index.time < time) &
                    (dHdl_restr_no_index.time >= time - dt)]
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
                    (dHdl_coul_no_index.time >= time - dt)]
                dHdl_local = dHdl_local.set_index(index_colnames)
                ti = TI().fit(dHdl_local)
                ti_coul_list.append(-self.conv_fac * ti.delta_f_.iloc[0, -1])
                ti_coul_sd_list.append(
                    self.conv_fac * ti.d_delta_f_.iloc[0, -1])

            # VDW
            if self.lambda_vdw:
                dHdl_local = dHdl_vdw_no_index[
                    (dHdl_vdw_no_index.time < time) &
                    (dHdl_vdw_no_index.time >= time - dt)]
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

        time_list = np.arange(dt, self.prod_time + dt, dt)
        self.convergence_data = {'time': time_list}

        for time in time_list:
            local_table = FreeEner.get_bar(
                self.xvg_file_list, begin_time=time - dt, end_time=time)
            tot_contrib = local_table['DG (kT)'].sum() * self.conv_fac
            tot_list.append(-tot_contrib)
            tot_std_contrib = sum(self.table['+/-']**2)**0.5 * self.conv_fac
            tot_sd_list.append(tot_std_contrib)

            if self.lambda_restr:
                restr_contrib = local_table['DG (kT)'][:restr_num - 1].sum()
                restr_contrib *= self.conv_fac
                restr_std_contrib = sum(
                    local_table['+/-'][:restr_num - 1]**2)**0.5
                restr_std_contrib *= self.conv_fac
                restr_list.append(-restr_contrib)
                restr_sd_list.append(restr_std_contrib)

            if self.lambda_coul:
                coul_contrib = local_table['DG (kT)'][
                    max(restr_num - 1, 0):restr_num + coul_num - 1].sum()
                coul_contrib *= self.conv_fac
                coul_std_contrib = sum(local_table['+/-'][
                    max(restr_num - 1, 0):restr_num + coul_num - 1]**2)**0.5
                coul_std_contrib *= self.conv_fac
                coul_list.append(-coul_contrib)
                coul_sd_list.append(coul_std_contrib)

            if self.lambda_vdw:
                vdw_contrib = local_table['DG (kT)'][
                    restr_num + coul_num - 1:].sum()
                vdw_contrib *= self.conv_fac
                vdw_std_contrib = sum(
                    local_table['+/-'][restr_num + coul_num - 1:]**2)**0.5
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
        width_ratio = 0.85

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
        axs[0].grid()
        box = axs[0].get_position()
        axs[0].set_position([box.x0, box.y0,
                             box.width * width_ratio,
                             box.height])

        ti_tot_list = self.convergence_data['tot'][0]
        ti_tot_sd_list = self.convergence_data['tot'][1]
        axs[1].errorbar(
            time_list, ti_tot_list, ti_tot_sd_list, label=r'$\Delta G_{tot}$')
        axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        axs[1].set_ylabel(r'$\Delta G \; ({})$'.format(self.unit_graph))
        axs[1].yaxis.set_label_coords(-0.1, 1)
        axs[1].set_xlabel(r'$time \; (ps)$')
        axs[1].grid()
        box = axs[1].get_position()
        axs[1].set_position([box.x0, box.y0,
                             box.width * width_ratio,
                             box.height])

        if graph_out:
            plt.savefig(graph_out)

        return(fig)

    @staticmethod
    def symmetry_correction(smile, temp=300):
        r""" Compute symmetry correction
        :math:`\Delta_{sym} = −k T ln(\sigma)`

        return value in kcal/mol

        .. code-block:: python

            > FreeEner.symmetry_correction('c1ccccc1')
            -1.4814...
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
                    f' sigma is set to {sigma}')

        GK = - KB * 1e-3 / 4.184  # Gas constant kcal/mol/K

        return GK * temp * math.log(sigma)

    """
    Test with threads:

    def run_threads(self, lambda_coul_list, lambda_vdw_list,
                    lambda_restr_list=[],
                    mbar=False, dir_name='free_ener_run',
                    em_steps=5000, nvt_time=10, npt_time=10, prod_time=100,
                    dt=0.002, temp=300.0, thread_num=4,
                    temp_groups='Protein non-Protein', maxwarn=1,
                    monitor_tool=monitor.PROGRESS_BAR):
        \"""Compute free energy to transfer a molecule from the
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

        \"""

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
                               'lincs_order': 8,  # Try to fix the LINCS Errors
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

        from concurrent.futures import ThreadPoolExecutor

        lambda_list = list(
            range(lambda_restr_num + lambda_coul_num + lambda_vdw_num))
        start_dir = os.path.abspath(".")

        def run_lambda(i):
            os.chdir(start_dir)
            print('change dir to:', start_dir)
            logger.info('Compute lambda {} / {}'.format(
                i + 1, lambda_restr_num + lambda_coul_num + lambda_vdw_num))

            lambda_sys = FreeEner.compute_lambda_point(
                self.gmxsys, i,
                self.mol_name,
                os.path.join(self.out_folder,
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
            return(lambda_sys)

        with ThreadPoolExecutor(max_workers=thread_num) as executor:
            lambda_results = executor.map(run_lambda, lambda_list)

        for result in lambda_results:
            self.lambda_sys_list.append(result)
            self.xvg_file_list.append(result.xvg)

        self.temp = temp
        self.lambda_coul = list(lambda_coul_list)
        self.lambda_vdw = list(lambda_vdw_list)
        self.lambda_restr = list(lambda_restr_list)
        self.prod_time = prod_time



    def extend_lambda_prod_threads(self, prod_time, thread_num=4):
        \"""
        \"""

        dt = float(self.lambda_sys_list[0].get_mdp_dict()['dt'])
        nsteps = prod_time / dt
        lambda_num = len(self.lambda_restr) + len(self.lambda_coul) +\
            len(self.lambda_vdw)
        xvg_file_list = []
        lambda_list = list(range(lambda_num))
        start_dir = os.path.abspath(".")

        def extend_lambda(lambda_i):
            os.chdir(start_dir)
            logger.info(f'Compute lambda {lambda_i+1} / {lambda_num}')

            self.lambda_sys_list[lambda_i].extend_sim(nsteps=nsteps)
            output = self.lambda_sys_list[lambda_i].get_all_output()
            return(output['xvg'])

        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=thread_num) as executor:
            lambda_results = executor.map(extend_lambda, lambda_list)

        for result in lambda_results:
            xvg_file_list += result

        self.prod_time = prod_time
        self.xvg_file_list = xvg_file_list
    """
