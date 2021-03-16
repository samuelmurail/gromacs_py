#!/usr/bin/env python3
# coding: utf-8
##################################
# #######   GROMACS 5   ##########
##################################


import os
import copy
import logging
import sys

from os_command_py import os_command
from pdb_manip_py import pdb_manip
from pdb_manip_py import pdb2pqr


from .topsys import TopSys
from .topmol import TopMol
from .itp import Itp

# Logging
logger = logging.getLogger(__name__)

# In case gmx5 is launched as main, relative import will failed
try:
    from ..tools import monitor, ambertools
except ImportError:
    logger.info("Relative import from .tools fails,"
                " use absolute import instead")
    import tools.monitor as monitor
    import tools.ambertools as ambertools


# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


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

    # Show log of top sys:
    from . import topsys
    topsys.show_log()

    if pdb_manip_log:
        # Show pdb_manip Logs:
        pdb_manip.show_log()


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
    if pdb_manip_log:
        # Show pdb_manip Logs:
        pdb_manip.show_log()


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
# Deduce the water gro file path
WATER_GRO = os.path.join(GMX_PATH, "share/gromacs/top/spc216.gro")

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

logger.info('FORCEFIELD_PATH_LIST = {}'.format(FORCEFIELD_PATH_LIST))

# Test folder path
GMX_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(GMX_LIB_DIR, "../test_files/")


# Global variable
# Protein heavy atoms
HA_NAME = ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'SG',
           'OG', 'OG1', 'CD', 'CD1', 'CD2', 'OD1', 'OD2', 'SD',
           'ND1', 'CE', 'CE1', 'CE2', 'CE3', 'OE1', 'OE2', 'NE',
           'NE1', 'NE2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'NH1',
           'NH2']
# DNA heavy atoms
HA_NAME += ['O5\'', 'C5\'', 'C4\'', 'O4\'', 'C1\'', 'N1', 'C6',
            'CG2', 'C5', 'C4', 'N4', 'N3', 'C2', 'O2', 'C3\'',
            'C2\'', 'O3\'', 'P', 'O1P', 'O2P', 'N9', 'C8', 'N7',
            'O6', 'N2', 'C7', 'N6', 'O4']

PROT_RES = ['GLY', 'HIS', 'HSP', 'HSE', 'HSD', 'HIP', 'HIE', 'HID',
            'ARG', 'LYS', 'ASP', 'ASPP', 'ASN', 'GLU', 'GLUP', 'GLN',
            'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', 'PRO', 'ALA',
            'ILE', 'PHE', 'TYR', 'TRP', 'VAL', 'LEU', 'MET']

DNARNA_RES = ['DA5', 'DA3', 'DAN', 'DA', 'DT5', 'DT3', 'DTN', 'DT', 'DC5',
              'DC3', 'DCN', 'DC', 'DG5', 'DG3', 'DGN', 'DG', 'RA5', 'RA3',
              'RAN', 'RA', 'RU5', 'RU3', 'RUN', 'RU', 'RC5', 'RC3', 'RCN',
              'RC', 'RG5', 'RG3', 'RGN', 'RG']

################################
# ## Gromacs System Object #####
################################


class GmxSys:
    """Gromacs system encapsulation class.

    This class can be used to launch most of gromacs
    commands (pdb2gmx, grompp, mdrun, trjconv, editconf, genconf, ...).
    After each steps, outputs file paths of gromacs commands are store in
    the class variable, like ``corr_file``, ``top_file``, tpr ...
    Most function will need the ``corr_file`` and/or ``top_file`` variable
    to be defined.

    The GmxSys object can be considered as a md simulation system. Each
    operation on the object will affect the object variables.

    The variables ``nt``, ``ntmpi`` and ``gpu_id`` are only used by functions
    which run simulations ``run_simulation()`` like ``em()`` or
    ``production()``.

    :param name: generic name of the system
    :type name: str

    :param sim_name: name of the simulation (used to create .tpr .log .edr ...)
    :type sim_name: str, optional

    :param coor_file: path of the coordinate file (.pdb, .gro)
    :type coor_file: str, optional

    :param top_file: path of the .top file
    :type top_file: str, optional

    :param tpr: path of the .tpr file
    :type tpr: str, optional

    :param mdp: path of the .mdp file
    :type mdp: str, optional

    :param xtc: path of the .xtc file
    :type xtc: str, optional

    :param edr: path of the .edr file
    :type edr: str, optional

    :param log: path of the .log file
    :type log: str, optional

    :param ndx: path of the .ndx file
    :type ndx: str, optional

    :param nt: Total number of threads to start
    :type nt: int, default=0

    :param ntmpi: Number of thread-MPI threads to start
    :type ntmpi: int, default=0

    :param gpu_id: List of GPU device id-s to use, specifies the per-node
        PP rank to GPU mapping
    :type gpu_id: str, default=None

    :param sys_history: List of previous GmxSys() states
    :type sys_history: list of GmxSys()

    :Example:

    .. code-block:: python

        > show_log()
        > TEST_OUT = str(getfixture('tmpdir'))
        > prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        > ###################################
        > ####   Create the topologie:   ###
        > ###################################
        > prot.prepare_top(out_folder=os.path.join(TEST_OUT, 'top_SH3'), \
    vsite='hydrogens') #doctest: +ELLIPSIS
        Succeed to read file .../test_files/1y0m.pdb ,  648 atoms found
        Succeed to read file .../test_files/1y0m.pdb ,  648 atoms found
        Succeed to save file tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method=propka \
    --with-ph=7.00 tmp_pdb2pqr.pdb 00_1y0m.pqr
        Succeed to read file 00_1y0m.pqr ,  996 atoms found
        Chain: A  Residue: 0 to 60
        Succeed to save file 01_1y0m_good_his.pdb
        - Create topologie
        gmx pdb2gmx -f 01_1y0m_good_his.pdb -o 1y0m_pdb2gmx.pdb -p \
    1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017 \
    -ignh -vsite hydrogens
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        > ###################################
        > ####    Add water and ions:     ###
        > ###################################
        > prot.solvate_add_ions(out_folder=os.path.join(TEST_OUT, 'top_sys'))\
        #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../top_SH3/1y0m_pdb2gmx.pdb -o \
    .../top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.1
        - Solvate the pbc box
        Copy topologie file and dependancies
        Copy topologie file and dependancies
        - Create the tpr file genion_1y0m_water_ion.tpr
        gmx grompp -f .../template/mini.mdp -c 1y0m_water.pdb -r \
    1y0m_water.pdb -p 1y0m_water_ion.top -po out_mini.mdp -o \
    genion_1y0m_water_ion.tpr -maxwarn 1
        - Add ions to the system with an ionic concentration of 0.15 M , \
    sytem charge = 0.0 water num= 4775
        Add ions : NA : 13   CL : 13
        gmx genion -s genion_1y0m_water_ion.tpr -p 1y0m_water_ion.top -o \
    1y0m_water_ion.gro -np 13 -pname NA -nn 13 -nname CL
        > ###################################
        > ####    Minimize the system     ###
        > ###################################
        > prot.em(out_folder=os.path.join(TEST_OUT, 'em_SH3'), nsteps=10, \
    constraints='none')
        - Create the tpr file 1y0m.tpr
        gmx grompp -f 1y0m.mdp -c ../top_sys/1y0m_water_ion.gro -r \
    ../top_sys/1y0m_water_ion.gro -p ../top_sys/1y0m_water_ion.top -po \
    out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
        - Launch the simulation 1y0m.tpr
        gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        > ###################################
        > ####    Create a D peptide      ###
        > ###################################
        > pep = GmxSys(name='D')
        > pep.create_peptide(sequence='D', out_folder=os.path.join(TEST_OUT, \
    'top_D'), em_nsteps=10, equi_nsteps=0, vsite='hydrogens')
        -Make peptide: D
        residue name:X
        residue name:D
        Succeed to save file .../top_D/D.pdb
        - Create topologie
        gmx pdb2gmx -f ../D.pdb -o D_pdb2gmx.pdb -p D_pdb2gmx.top -i \
D_posre.itp -water tip3p -ff charmm36-jul2017 -ignh -ter -vsite hydrogens
        Molecule topologie present in D_pdb2gmx.top , extract the topologie \
    in a separate file: D_pdb2gmx.itp
        Protein_chain_P
        - ITP file: D_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_P
        Rewrite topologie: D_pdb2gmx.top
        - Create pbc box
        gmx editconf -f .../top_D/00_top/D_pdb2gmx.pdb -o \
    .../top_D/00_top/D_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        - Create the tpr file D.tpr
        gmx grompp -f D.mdp -c ../00_top/D_pdb2gmx_box.pdb -r \
    ../00_top/D_pdb2gmx_box.pdb -p ../00_top/D_pdb2gmx.top -po out_D.mdp -o \
    D.tpr -maxwarn 1
        - Launch the simulation D.tpr
        gmx mdrun -s D.tpr -deffnm D -nt 0 -ntmpi 0 -nsteps -2 -nocopyright
        > #######################################################
        > ### Insert 4 copy of the peptide in the SH3 system: ###
        > #######################################################
        > prot.insert_mol_sys(mol_gromacs=pep, mol_num=4, new_name='SH3_D', \
    out_folder=os.path.join(TEST_OUT, 'top_D_SH3')) #doctest: +ELLIPSIS
        - Convert trj/coor
        gmx trjconv -f ...D.gro -o ...D_compact.pdb -s ...D.gro -ur compact \
    -pbc none
        Succeed to read file ...D_compact.pdb ,  22 atoms found
        - Copy pbc box using genconf
        Succeed to read file ...D_compact_copy_box.pdb ,  88 atoms found
        Succeed to save file ...D_compact_copy_box.pdb
        Res num: 8
        - Convert trj/coor
        gmx trjconv -f ../em_SH3/1y0m.gro -o ../em_SH3/1y0m_compact.pdb -s \
    ../em_SH3/1y0m.tpr -ur compact -pbc mol
        Concat files: ['../em_SH3/1y0m_compact.pdb', \
    '../top_D/01_mini/D_compact_copy_box.pdb']
        Succeed to save concat file:  SH3_D_pre_mix.pdb
        Succeed to read file SH3_D_pre_mix.pdb ,  15425 atoms found
        Insert mol in system
        Residue list = [4836, 4837, 4838, 4839, 4840, 4841, 4842, 4843]
        Insert 4 mol of 2 residues each
        insert mol   1, water mol   ..., time=0...
        Warning atom 1MCH mass could not be founded
        Warning atom 2MCH mass could not be founded
        Warning atom 1HH3 mass could not be founded
        Warning atom 2HH3 mass could not be founded
        Warning atom 3HH3 mass could not be founded
        insert mol   2, water mol   ..., time=0...
        Warning atom 1MCH mass could not be founded
        Warning atom 2MCH mass could not be founded
        Warning atom 1HH3 mass could not be founded
        Warning atom 2HH3 mass could not be founded
        Warning atom 3HH3 mass could not be founded
        insert mol   3, water mol   ..., time=0...
        Warning atom 1MCH mass could not be founded
        Warning atom 2MCH mass could not be founded
        Warning atom 1HH3 mass could not be founded
        Warning atom 2HH3 mass could not be founded
        Warning atom 3HH3 mass could not be founded
        insert mol   4, water mol   ..., time=0...
        Warning atom 1MCH mass could not be founded
        Warning atom 2MCH mass could not be founded
        Warning atom 1HH3 mass could not be founded
        Warning atom 2HH3 mass could not be founded
        Warning atom 3HH3 mass could not be founded
        Delete ... overlapping water atoms
        Succeed to save file SH3_D.pdb
        Ligand
        Add 4 mol ...D_pdb2gmx.itp
        Succeed to read file SH3_D.pdb ,  15... atoms found
        Water num: 47...
        [{'name': 'Protein_chain_A', 'num': '1'}, {'name': 'SOL', ...\
    {'name': 'Ligand', 'num': '4'}]
        CHARGE: -4.0
        Should neutralize the system
        Copy topologie file and dependancies
        - Create the tpr file genion_SH3_D_neutral.tpr
        gmx grompp -f .../template/mini.mdp -c SH3_D.pdb -r SH3_D.pdb -p \
    SH3_D_neutral.top -po out_mini.mdp -o genion_SH3_D_neutral.tpr -maxwarn 1
        - Add ions to the system with an ionic concentration of 0 M , \
    sytem charge = -4.0 water num= 47...
        Add ions : NA : 4   CL : 0
        gmx genion -s genion_SH3_D_neutral.tpr -p SH3_D_neutral.top -o \
    SH3_D_neutral.gro -np 4 -pname NA -nn 0 -nname CL
        > ################################
        > ####   Minimize the system   ###
        > ################################
        > prot.em_2_steps(out_folder=os.path.join(TEST_OUT, 'top_D_SH3'), \
    no_constr_nsteps=10, constr_nsteps=10)
        - Create the tpr file Init_em_1y0m.tpr
        gmx grompp -f Init_em_1y0m.mdp -c SH3_D_neutral.gro -r \
SH3_D_neutral.gro -p SH3_D_neutral.top -po out_Init_em_1y0m.mdp \
-o Init_em_1y0m.tpr -maxwarn 1
        - Launch the simulation Init_em_1y0m.tpr
        gmx mdrun -s Init_em_1y0m.tpr -deffnm Init_em_1y0m -nt 0 -ntmpi 0 \
    -nsteps -2 -nocopyright
        - Create the tpr file 1y0m.tpr
        gmx grompp -f 1y0m.mdp -c Init_em_1y0m.gro -r Init_em_1y0m.gro -p \
    SH3_D_neutral.top -po out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
        - Launch the simulation 1y0m.tpr
        gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        > ##################################
        > ####    Show system history    ###
        > ##################################
        > prot.display_history() #doctest: +ELLIPSIS
        State -3:
        <BLANKLINE>
        name         : 1y0m
        sim_name     : genion_1y0m_water_ion
        coor_file    : .../top_sys/1y0m_water_ion.gro
        top_file     : .../top_sys/1y0m_water_ion.top
        tpr          : .../top_sys/genion_1y0m_water_ion.tpr
        mdp          : ...template/mini.mdp
        nt           : 0
        ntmpi        : 0
        sys_history  : 0
        <BLANKLINE>
        State -2:
        <BLANKLINE>
        name         : 1y0m
        sim_name     : genion_SH3_D_neutral
        coor_file    : .../top_D_SH3/SH3_D_neutral.gro
        top_file     : .../top_D_SH3/SH3_D_neutral.top
        tpr          : .../top_D_SH3/genion_SH3_D_neutral.tpr
        mdp          : ...template/mini.mdp
        xtc          : .../em_SH3/1y0m.trr
        edr          : .../em_SH3/1y0m.edr
        log          : .../em_SH3/1y0m.log
        nt           : 0
        ntmpi        : 0
        sys_history  : 0
        <BLANKLINE>
        State -1:
        <BLANKLINE>
        name         : 1y0m
        sim_name     : Init_em_1y0m
        coor_file    : .../top_D_SH3/Init_em_1y0m.gro
        top_file     : .../top_D_SH3/SH3_D_neutral.top
        tpr          : .../top_D_SH3/Init_em_1y0m.tpr
        mdp          : .../top_D_SH3/Init_em_1y0m.mdp
        xtc          : .../top_D_SH3/Init_em_1y0m.trr
        edr          : .../top_D_SH3/Init_em_1y0m.edr
        log          : .../top_D_SH3/Init_em_1y0m.log
        nt           : 0
        ntmpi        : 0
        sys_history  : 0
        <BLANKLINE>
        > ###################################
        > ####   Equilibrate the system   ###
        > ###################################
        > equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME, \
    "template/equi_vsites.mdp")
        > mdp_options = {'nsteps': 100, 'define': '-DPOSRES', 'dt': 0.001}
        > prot.run_md_sim(out_folder=os.path.join(TEST_OUT, 'equi_HA_D_SH3'), \
    name="equi_HA_D_SH3", mdp_template=equi_template_mdp,\
                            mdp_options=mdp_options)
        - Create the tpr file equi_HA_D_SH3.tpr
        gmx grompp -f equi_HA_D_SH3.mdp -c ../top_D_SH3/1y0m.gro -r \
    ../top_D_SH3/1y0m.gro -p ../top_D_SH3/SH3_D_neutral.top -po \
    out_equi_HA_D_SH3.mdp -o equi_HA_D_SH3.tpr -maxwarn 0
        - Launch the simulation equi_HA_D_SH3.tpr
        gmx mdrun -s equi_HA_D_SH3.tpr -deffnm equi_HA_D_SH3 -nt 0 -ntmpi 0 \
    -nsteps -2 -nocopyright
        > prot.get_simulation_time() #doctest: +ELLIPSIS
        - Get simulation time from : .../equi_HA_D_SH3/equi_HA_D_SH3.cpt
        gmx check -f .../equi_HA_D_SH3/equi_HA_D_SH3.cpt
        0.1
        > prot.convert_trj(traj=False) #doctest: +ELLIPSIS
        - Convert trj/coor
        gmx trjconv -f .../equi_HA_D_SH3/equi_HA_D_SH3.gro -o \
    .../equi_HA_D_SH3/equi_HA_D_SH3_compact.pdb -s \
    .../equi_HA_D_SH3/equi_HA_D_SH3.tpr -ur compact -pbc mol
        > prot.display() #doctest: +ELLIPSIS
        name         : 1y0m
        sim_name     : equi_HA_D_SH3
        coor_file    : .../equi_HA_D_SH3/equi_HA_D_SH3_compact.pdb
        top_file     : .../top_D_SH3/SH3_D_neutral.top
        tpr          : .../equi_HA_D_SH3/equi_HA_D_SH3.tpr
        mdp          : .../equi_HA_D_SH3/equi_HA_D_SH3.mdp
        xtc          : .../equi_HA_D_SH3/equi_HA_D_SH3.xtc
        edr          : .../equi_HA_D_SH3/equi_HA_D_SH3.edr
        log          : .../equi_HA_D_SH3/equi_HA_D_SH3.log
        nt           : 0
        ntmpi        : 0
        sys_history  : 4
        > #########################################
        > ### Extract Potential Energy and Temp ###
        > #########################################
        > ener_pd = prot.get_ener(['Potential', 'Temp'])  #doctest: +ELLIPSIS
        - Extract energy
        gmx energy -f .../equi_HA_D_SH3/equi_HA_D_SH3.edr -o tmp_edr.xvg
        > ener_pd['Potential'].mean() #doctest: +ELLIPSIS
        -2...
        > rmsd_pd = prot.get_rmsd(['C-alpha', 'Protein'])  #doctest: +ELLIPSIS
        - Extract RMSD
        - Create the ndx file ...equi_HA_D_SH3.ndx
        gmx make_ndx -f ...equi_HA_D_SH3_compact.pdb -o ...equi_HA_D_SH3.ndx
        gmx rms -s ...equi_HA_D_SH3.tpr -f ...equi_HA_D_SH3.xtc -n ...\
    equi_HA_D_SH3.ndx -o tmp_rmsd.xvg -fit rot+trans -ng 1 -pbc no
        > rmsd_pd #doctest: +ELLIPSIS
           time ...Protein
        0   0.0...
        > rmsf_pd = prot.get_rmsf(['Protein'], res="yes")  #doctest: +ELLIPSIS
        - Extract RMSF
        gmx rmsf -s ...equi_HA_D_SH3.tpr -f ...equi_HA_D_SH3.xtc -n ...\
    equi_HA_D_SH3.ndx -o tmp_rmsf.xvg -fit no -res yes
        > rmsf_pd #doctest: +ELLIPSIS
            Residue    RMSF
        0       791  ...

    .. note::
        An history of all command used could be saved.

    """

    def __init__(self, name=None, coor_file=None, top_file=None, tpr=None):
        """**__init__:**
        All files path are check when the variable are assigned with the
        command ``os_command.full_path_and_check()``, if the file does not
        exist, the program will crash. This is usefull as when gromacs
        commands fail, the only way to realise it, is to check the created
        files.

        :param name: generic name of the system
        :type name: str, optional

        :param coor_file: path of the coor file
        :type coor_file: str, optional

        :param tpr: path of the tpr file
        :type tpr: str, optional

        .. note::
            Not sure that coor_file, top_file and tpr should the only
            parameters.

        .. todo::
            Add exeption when ``os_command.full_path_and_check()`` crash
            because files do not exist.
        """

        # System name:
        self.name = name
        self.sim_name = None

        # Simulation inputs:
        self.coor_file = coor_file
        self.top_file = top_file
        self.tpr = tpr
        self._ndx = None

        # Simulation outputs:
        self._mdp = None
        self._xtc = None
        self._edr = None
        self._log = None
        self._xvg = None

        # mdrun default run values:
        self.nt = 0
        self.ntmpi = 0
        self.gpu_id = None

        # System information history:
        # list of GmxSys()
        self.sys_history = []

    @property
    def coor_file(self):
        if self._coor_file is not None:
            return os.path.relpath(self._coor_file)
        return None

    @property
    def top_file(self):
        if self._top_file is not None:
            return os.path.relpath(self._top_file)
        return None

    @property
    def tpr(self):
        if self._tpr is not None:
            return os.path.relpath(self._tpr)
        return None

    @property
    def mdp(self):
        if self._mdp is not None:
            return os.path.relpath(self._mdp)
        return None

    @property
    def xtc(self):
        if self._xtc is not None:
            return os.path.relpath(self._xtc)
        return None

    @property
    def edr(self):
        if self._edr is not None:
            return os.path.relpath(self._edr)
        return None

    @property
    def log(self):
        if self._log is not None:
            return os.path.relpath(self._log)
        return None

    @property
    def xvg(self):
        if self._xvg is not None:
            return os.path.relpath(self._xvg)
        return None

    @property
    def ndx(self):
        if self._ndx is not None:
            return os.path.relpath(self._ndx)
        return None

    @coor_file.setter
    def coor_file(self, coor_file):
        if coor_file is not None:
            self._coor_file = os_command.full_path_and_check(coor_file)
        else:
            self._coor_file = None

    @top_file.setter
    def top_file(self, top_file):
        if top_file is not None:
            self._top_file = os_command.full_path_and_check(top_file)
        else:
            self._top_file = None

    @tpr.setter
    def tpr(self, tpr):
        if tpr is not None:
            self._tpr = os_command.full_path_and_check(tpr)
        else:
            self._tpr = None

    @mdp.setter
    def mdp(self, mdp):
        self._mdp = os_command.full_path_and_check(mdp)

    @xtc.setter
    def xtc(self, xtc):
        self._xtc = os_command.full_path_and_check(xtc)

    @edr.setter
    def edr(self, edr):
        self._edr = os_command.full_path_and_check(edr)

    @log.setter
    def log(self, log):
        self._log = os_command.full_path_and_check(log)

    @xvg.setter
    def xvg(self, xvg):
        self._xvg = os_command.full_path_and_check(xvg)

    @ndx.setter
    def ndx(self, ndx):
        self._ndx = os_command.full_path_and_check(ndx)

    def display(self):
        """Display defined attribute of the GmxSys object.
        """
        # print("Coor : ", self.coor_file, "\nTop : ", self.top_file)

        # Order dict is only necessary for python 3.5, where dict are not
        # ordered.
        # This is only require for doctest
        numbermap = {'name': 1,
                     'sim_name': 2,
                     '_coor_file': 3,
                     '_top_file': 4,
                     '_tpr': 5,
                     '_ndx': 6,
                     '_mdp': 7,
                     '_xtc': 8,
                     '_edr': 9,
                     '_log': 10,
                     '_xvg': 11,
                     'nt': 12,
                     'ntmpi': 13,
                     'gpu_id': 14,
                     'sys_history': 15}

        attr_list = [attr for attr in vars(self) if not attr.startswith('__')]
        for attr in sorted(attr_list, key=numbermap.__getitem__):
            if attr[0] == "_":
                to_show = attr[1:]
            else:
                to_show = attr
            if attr == 'sys_history':
                print("{:12} : {}".format(
                    to_show, len(getattr(self, to_show))))
            elif getattr(self, to_show) is not None:
                print("{:12} : {}".format(to_show, getattr(self, to_show)))

    def save_state(self):
        """ Save last state
        """

        prev_state = copy.deepcopy(self)
        prev_state.sys_history = []
        self.sys_history.append(prev_state)

    def display_history(self):
        """ Show all history
        """

        for i, history in enumerate(self.sys_history):
            print("State {}:\n".format(i - len(self.sys_history)))
            history.display()
            print()

    def view_coor(self):
        """ Return a `nglview` object to view the object coordinates
        in a jupyter notebook with the module ``nglview``.

        :Example:

        >>> import nglview as nv #doctest: +SKIP
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> view = prot.view_coor() #doctest: +SKIP
        >>> view #doctest: +SKIP
        """
        try:
            import nglview as nv
        except ImportError:
            logger.warning(
                'Could not load nglview \nInstall it using conda:\n' +
                'conda install -c conda-forge nglview')
            return

        if self.coor_file[-3:] not in ['pdb', 'gro']:
            if self.tpr is None:
                self.tpr = self.coor_file
            self.convert_trj(traj=False, pbc='none')

        coor = pdb_manip.Coor(self.coor_file)

        struct_str = nv.TextStructure(coor.get_structure_string())
        return nv.NGLWidget(struct_str)

    def view_traj(self):
        """ Return a `nglview` object to view the object trajectorie
        in a jupyter notebook with the module ``nglview``.

        :Example:

        >>> import nglview as nv #doctest: +SKIP
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> view = prot.view_traj() #doctest: +SKIP
        >>> view #doctest: +SKIP

        .. note::
            This function has a dependencies with the ``simpletraj`` a
            lightweight coordinate-only trajectory reader. Install it
            using ``pip install simpletraj`` or ``conda install simpletraj``.
        """

        import nglview as nv

        simple_traj = nv.SimpletrajTrajectory(self.xtc, self.coor_file)
        return nv.NGLWidget(simple_traj)

    #########################################################
    # ###########  TOPOLOGIE RELATED FUNCTIONS  #############
    #########################################################

    def add_top(self, out_folder, name=None, ff="charmm36-jul2017",
                water="tip3p", check_file_out=True, pdb2gmx_option_dict={},
                input_pdb2gmx="", posre_post=""):
        """Launch the pdb2gmx command.

        The objet variable self.coor_file has to be defined before launching
        this function. ``pdb2gmx`` will create a new coordinate file
        ``name+"_pdb2gmx.pdb"``, a topologie ``name+"_pdb2gmx.top"`` and
        several molecule itp and posre files.
        If name is not defined, it will use the object name.

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ff: forcefield
        :type ff: str, optional, default="charmm36"

        :param water: water model
        :type water: str, optional, default="tip3p"

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param pdb2gmx_option_dict: dictionnary of option for pdb2gmx, for
            example if you want to ignore input hydrogens use:
            ``{'ignh': None}``. The '-' before the option is to avoid.
        :type pdb2gmx_option_dict: dict, optional, default=None

        :param input_pdb2gmx: input for pdb2gmx request
        :type input_pdb2gmx: str, optional, default=None

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> #Basic usage :
        >>> prot.add_top(out_folder=TEST_OUT+'/add_top/top_SH3')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> #########################################
        >>> # Use of different options for pdb2gmx: #
        >>> #########################################
        >>> # Ignore hydrogens: 'ignh': None
        >>> # Define amino acid termini: 'ter': None
        >>> # Needs to answer pdb2gmx request concerning termini
        >>> # with: input_pdb2gmx ="1 \\n 0"
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/add_top/top_SH3_2/',
        ...     pdb2gmx_option_dict={'ignh': None, 'ter': None},
        ...     input_pdb2gmx="1 \\n 0") #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017 \
-ignh -ter
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top

        .. note::
            To avoid conflict with focefields, the environment variable
            ``$GMXLIB`` is change to ``GROMACS_MOD_DIRNAME+"/template/"``
            where curently only charmm36 is present, if you want to use
            another forcefield, copy your forcefield folder in
            ``GROMACS_MOD_DIRNAME+"/template/"``, or change the current code.
        """

        logger.info("- Create topologie")

        # Get absolute path:
        if self.coor_file is None:
            raise RuntimeError('self.coor_file not defined')

        start_dir = os.path.abspath(".")

        # If name is not define use the object name
        if name is None:
            name = self.name

        # Create and go in out_folder:
        # This is necessary for the topologie creation
        os_command.create_and_go_dir(out_folder)

        new_coor = name + "_pdb2gmx.pdb"
        top_file = name + "_pdb2gmx.top"
        posre_file = name + "_posre.itp"

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(new_coor):
            logger.info("create_top not launched {} already exist".format(
                new_coor))
            self.coor_file = new_coor
            self.top_file = top_file
            os.chdir(start_dir)
            return

        # Define pdb2gmx command:
        cmd_top = os_command.Command([GMX_BIN, "pdb2gmx",
                                      "-f", self.coor_file,
                                      "-o", new_coor,
                                      "-p", top_file,
                                      "-i", posre_file,
                                      "-water", water,
                                      "-ff", ff], **pdb2gmx_option_dict)

        # Define the forcefield dir to the env GMXLIB,
        cmd_top.define_env(my_env=os.environ.update({'GMXLIB': GMXLIB_var}))
        cmd_top.display()
        cmd_top.run(input_pdb2gmx)

        # First read and save to fix the molecule top include in .top:
        top = TopSys(top_file)
        top.write_file(top_file)

        # Now it can add posre files properly:
        top = TopSys(top_file)
        top.add_posre(posre_name="POSRES_HA_LOW" + posre_post, selec_dict={
            'atom_name': HA_NAME},
            fc=[100, 100, 100])
        top.add_posre(posre_name="POSRES_CA_LOW" + posre_post, selec_dict={
            'atom_name': ['CA', 'P']},
            fc=[100, 100, 100])
        top.add_posre(posre_name="POSRES_CA" + posre_post, selec_dict={
            'atom_name': ['CA', 'P']},
            fc=[1000, 1000, 1000])

        self.coor_file = new_coor
        self.top_file = top_file

        os.chdir(start_dir)

    def prepare_top(self, out_folder, name=None,
                    vsite="none", ignore_ZN=True,
                    ff="charmm36-jul2017", ph=7.0,
                    res_prot_dict=None,
                    include_mol={}):
        """Prepare the topologie of a protein:

        1. compute hisdine protonation with ``pdb2pqr``.

        2. Change Histidine resname according to the protonation.

        3. Correct cystein resname.

        4. Correct chain ID's.

        5. Zinc Finger: Add Zinc in the pdb and change residue type of
            CYS and HIS coordinating the Zinc.

        6. Finally compute the topologie with pdb2gmx add_top().

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param vsite: option for topologie's bonds constraints ("none",
            "hydrogens", "aromatics")
        :type vsite: str, optional, default="none"

        :param ignore_ZN: option for not adding parameters to ZINC finger
        :type ignore_ZN: bool, optional, default=False

        :param ff: forcefield
        :type ff: str, optional, default="charmm36-jul2017"

        :param ph: pH to assign AA protonation (using pdb2pqr)
        :type ph: float, optional, default=7.0

        :param res_prot_dict: option to define manually protonation
        :type res_prot_dict: dict, optional, default=None

        :param include_mol: list of ligand's residue name to include
        :type include_mol: list, optional, default=[]

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> # Create the topologie of a protein and do a minimisation:
        >>> dna_lig = GmxSys(name='1D30', coor_file=TEST_PATH+'/1D30.pdb')
        >>> dna_lig.prepare_top(out_folder=TEST_OUT+'/prepare_top/top_dna/', \
ff='amber99sb-ildn', include_mol={'DAP':\
'NC(=N)c1ccc(cc1)c2[nH]c3cc(ccc3c2)C(N)=N'}) #doctest: +ELLIPSIS
        Succeed to read file .../test_files/1D30.pdb ,  532 atoms found
        Succeed to read file .../test_files/1D30.pdb ,  532 atoms found
        Succeed to save file tmp_pdb2pqr.pdb
        pdb2pqr... --ff AMBER --ffout AMBER --chain --ph-calc-method=propka \
--with-ph=7.00 tmp_pdb2pqr.pdb 00_1D30.pqr
        Succeed to read file 00_1D30.pqr ,  758 atoms found
        Succeed to read file .../test_files/1D30.pdb ,  532 atoms found
        Succeed to save file DAP.pdb
        Succeed to read file DAP.pdb ,  21 atoms found
        Succeed to save file DAP_0.pdb
        Succeed to read file DAP_0.pdb ,  21 atoms found
        Succeed to read file DAP_0_h.pdb ,  36 atoms found
        Succeed to save file DAP_0_h.pdb
        Succeed to read file DAP_0_h.pdb ,  36 atoms found
        Succeed to save file DAP_0_h.pdb
        Succeed to save concat file:  DAP_h.pdb
        Succeed to read file DAP_h.pdb ,  36 atoms found
        Succeed to save file DAP_h_unique.pdb
        acpype... -i DAP_h_unique.pdb -b DAP -c bcc -a gaff -o gmx -n 0
        DAP
        Succeed to save file 01_1D30_good_his.pdb
        - Create topologie
        gmx pdb2gmx -f 01_1D30_good_his.pdb -o 1D30_pdb2gmx.pdb -p \
1D30_pdb2gmx.top -i 1D30_posre.itp -water tip3p -ff amber99sb-ildn \
-ignh -vsite none
        Add Molecule: DAP
        Add 1 mol .../top_dna/DAP.acpype/DAP_GMX.itp
        Concat files: ['1D30_pdb2gmx.pdb', '.../top_dna/DAP_h.pdb']
        Succeed to save concat file:  1D30_pdb2gmx_mol.pdb
        >>> dna_lig.em(out_folder=TEST_OUT + '/prepare_top/em_dna', \
nsteps=10, create_box_flag=True, maxwarn=1) #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../top_dna/1D30_pdb2gmx_mol.pdb -o \
.../top_dna/1D30_pdb2gmx_mol_box.pdb -bt dodecahedron -d 1.0
        - Create the tpr file 1D30.tpr
        gmx grompp -f 1D30.mdp -c .../top_dna/1D30_pdb2gmx_mol_box.pdb -r \
.../top_dna/1D30_pdb2gmx_mol_box.pdb -p .../top_dna/1D30_pdb2gmx_mol.top -po \
out_1D30.mdp -o 1D30.tpr -maxwarn 1
        - Launch the simulation 1D30.tpr
        gmx mdrun -s 1D30.tpr -deffnm 1D30 -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        >>> lig = dna_lig.extract_mol_sys(out_folder=TEST_OUT+\
'/prepare_top/top_DAP/', res_name='DAP') #doctest: +ELLIPSIS
        - Convert trj/coor
        gmx trjconv -f ...1D30.gro -o ...1D30_compact.pdb -s ...1D30.tpr \
-ur compact -pbc mol
        Succeed to read file ...1D30_compact.pdb ,  794 atoms found
        Succeed to save file ...DAP_only.pdb
        Forcefield include :
         amber99sb-ildn
        - ITP file: DAP_GMX_atomtypes
        - molecules defined in the itp file:
        - ITP file: tip3p
        - molecules defined in the itp file:
        * SOL
        - ITP file: DAP_GMX
        - molecules defined in the itp file:
        * DAP
        Mol List:
           * 1 DAP
        Mol Name:
         Protein
        >>> lig.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f ...DAP_only.pdb -o ...DAP_only_box.pdb -bt \
dodecahedron -d 1.0
        >>> lig.solvate_box(out_folder=TEST_OUT + '/prepare_top/water_lig_top')
        - Solvate the pbc box
        Copy topologie file and dependancies
        >>> lig.em(out_folder=TEST_OUT + '/prepare_top/em_water_DAP', \
nsteps=10, maxwarn=1) #doctest: +ELLIPSIS
        - Create the tpr file DAP.tpr
        gmx grompp -f DAP.mdp -c ...DAP_water.pdb -r ...DAP_water.pdb -p \
...DAP_water.top -po out_DAP.mdp -o DAP.tpr -maxwarn 1
        - Launch the simulation DAP.tpr
        gmx mdrun -s DAP.tpr -deffnm DAP -nt 0 -ntmpi 0 -nsteps -2 -nocopyright

        .. note::
            No options are allowed (forcefield, water model, termini capping)
            except for vsites.

        """

        start_dir = os.path.abspath(".")

        # Create and go in out_folder:
        # This is necessary for the topologie creation
        os_command.create_and_go_dir(out_folder)

        # If name is not define use the object name
        if name is None:
            name = self.name

        # Save initial pdb file:
        start_pdb = self.coor_file
        start_coor = pdb_manip.Coor(start_pdb)

        # Get initial list of residue name
        res_name_input = start_coor.get_attribute_selection(
            attribute='res_name')

        if ff.startswith('amber'):
            pdb2pqr_ff = 'AMBER'
        else:
            pdb2pqr_ff = 'CHARMM'

        # Compute protonation:
        pdb2pqr.compute_pdb2pqr(self.coor_file,
                                "00_" + name + ".pqr",
                                ff=pdb2pqr_ff, ph=ph,
                                check_file_out=True)

        # Correct His resname
        coor_in = pdb_manip.Coor("00_" + name + ".pqr")
        res_name_pdb2pqr = coor_in.get_attribute_selection(
            attribute='res_name')
        res_name_pdb2pqr += ['HIS', 'HOH']

        mol_sys_list = []
        for resname in res_name_input:
            if resname not in res_name_pdb2pqr:
                if resname in include_mol:
                    # mol_top = ambertools.make_amber_top_mol(
                    #        start_pdb, resname, charge=include_mol[resname])
                    mol_top = ambertools.make_amber_top_mol_rdkit(
                            start_pdb, res_name=resname,
                            smile=include_mol[resname],
                            remove_h=True)
                    mol_top['name'] = resname
                    mol_sys_list.append(mol_top)
                else:
                    logger.warning("residue(s) {} not included,"
                                   " add this residue in ...".format(
                                    resname))

        coor_in.correct_water_name()
        if start_coor.get_aa_num() > 0:
            coor_in.correct_his_name()
            coor_in.correct_cys_name()
            coor_in.correct_chain()
            coor_in.correct_protonated_res()
            if res_prot_dict is not None:
                GmxSys.set_coor_aa_prot(coor_in, res_prot_dict,
                                        ff=pdb2pqr_ff)

        # NOTE The zinc finger should be removed from here:
        if not ignore_ZN:
            coor_in.add_zinc_finger(start_pdb)
        # else:
        #    zinc_in = False

        coor_in.write_pdb(pdb_out="01_" + name + "_good_his.pdb")

        self.coor_file = "01_" + name + "_good_his.pdb"

        # Compute topology options
        pdb2gmx_option_dict = {'vsite': vsite, 'ignh': None}
        # Options for system with zinc

        # if zinc_in:  # No more necessay if using octahedral dummy Zinc model
        #    pdb2gmx_option_dict['merge'] = 'all'

        self.add_top(out_folder=".", check_file_out=True, ff=ff,
                     pdb2gmx_option_dict=pdb2gmx_option_dict)

        # Add molecules coordinates and topologies:

        # Get the system topologie:
        sys_topologie = TopSys(self.top_file)

        pdb_mol_list = []
        # Add the molecule in the sys topologie and update the water num:
        for mol in mol_sys_list:
            logger.info('Add Molecule: {}'.format(mol['name']))

            pdb_mol_list.append(mol['coor'])
            # Add topologie:
            # mol['GmxSys'].display()
            mol_top = TopSys(mol['GmxSys'].top_file)
            mol_itp = mol_top.get_include_no_posre_file_list()
            sys_topologie.add_mol(mol_name=mol['name'],
                                  mol_itp_file=mol_itp[-1],
                                  mol_num=mol['num'])
            # Add atomtypes itp:
            sys_topologie.add_atomtypes(mol_itp[0])

        if mol_sys_list:
            # Add coordinates:
            GmxSys.concat_coor(self.coor_file, *pdb_mol_list,
                               pdb_out=self.coor_file[:-4] + '_mol.pdb',
                               check_file_out=False)
            self.coor_file = self.coor_file[:-4] + '_mol.pdb'
            # Save topologie
            sys_topologie.write_file(self.top_file[:-4] + '_mol.top')
            self.top_file = self.top_file[:-4] + '_mol.top'

        os.chdir(start_dir)

    def create_itp_atomtype_ion_octa_dummy(self, atomtypes,
                                           ion_name=['MN', 'ZN']):
        r"""
        Forcefield A and B values taken from :
        Duarte et al. 2014 J. Phys. Chem. B

        https://en.wikipedia.org/wiki/Lennard-Jones_potential

        :math:`A = 4 \epsilon \sigma^{12}`

        :math:`B = 4 \epsilon \sigma^6`

        :math:`\sigma = \sqrt[6]{\frac{A}{B}}`

        :math:`\epsilon = \frac{B^2}{4A}`

        ion_name=['NI', 'CO', 'ZN', 'MN', 'FE', 'MG', 'CA']

        ; MM  171         35
        ; D   0.05         0
        ;
        ; MN
        ; C12, C6 = 171**2, 35**2
        ; sigma = (C12/C6)**(1/6) = 1.697 A = 0.1697 nm
        ; eps = C6**2/(4*C12) = 12.829 Kcal mol-1 = 418.4 KJ mol -1
        ;
        ; DMN
        ; C12, C6 = 0, 35**2

        ; MN  25    36.938000    0.000  A   0.1697    12.829
        ; Gromacs unit is : kJ mol−1 nm−2   kb = kb_cal_A * 4.184 * 100

        ; Aqvist and Warshel JACS 1990
        ;
        ; MM  145        25
        ; D   0           0
        ; Kb = 1600 (kcal mol−1Å−2) and Kθ = 250 (kcal mol−1rad−2) and
        ; no bond between dummies.
        ; C12, C6 = 145**2, 25**2
        ; sigma = (C12/C6)**(1/6) = 1.7967 A = 0.17967 nm
        ; eps = C6**2/(4*C12) = 4.645 Kcal mol-1 = 19.4337 KJ mol -1
        ; MN  25    36.938000    0.000  A   0.17967    19.4337

        """
        # Add Topologie:

        # Compute sigma and epsilon form A and B
        for name, atom in atomtypes.items():
            C6, C12 = atom['A']**2, atom['B']**2
            # A in kcal.mol-1.Å(-12)
            # Convert to kJ.mol-1.nm(-12)
            C6 *= 4.184 * 10**-12
            # B in kcal.mol.Å(-3)
            # Convert to kJ.mol-1.nm(-6)
            C12 *= 4.184 * 10**-6
            if C12 != 0.0:
                sigma = (C6 / C12)**(1 / 6)
                eps = C12**2 / (4 * C6)
            else:
                # Is it a good fix ?, need to adjust sigma
                # to have high enough epsilon
                sigma = -0.1
                eps = C6 / (4 * sigma**(12))

            atom['eps'] = eps  # nm
            atom['sigma'] = sigma  # KJ.mol^(-1)

            # print(f'C6 ={(C6 * 4.184 * 1e-6):8.2e} '
            #       'C12 = {(C12 * 4.184 * 1e-12):7.4e}')
            print(f'{name:4} sig = {(sigma):7.4f} nm'
                  f' eps ={eps:8.2f} KJ.mol-1.nm-2')

        # Atomtype itp file:

        ion_atom_type = self.top_file[:-4] + '_atomtypes.itp'

        # Write the atomtype in a separate file:
        with open(ion_atom_type, 'w') as file:
            file.write('[ atomtypes ]\n')

            for name in ion_name+['D']:
                atom_dict = atomtypes[name]
                file.write(' {:3}      {:3}         {:8.5f}  {:8.5f}'
                           '   {}     {:12.5e}   {:12.5e}\n'.format(
                            name,
                            atom_dict['atnum'],
                            atom_dict['mass'],
                            atom_dict['charge'],
                            'A',
                            atom_dict['sigma'],
                            atom_dict['eps']))
        return ion_atom_type

    def create_itp_ion_octa_dummy(self, atomtypes, ion_name=['MN', 'ZN']):
        """
        """

        # Unit: distance in (nm), k in kcal.mol-1.Å-2
        bond_type = {'ID': {'d': 0.0900, 'k': 800},
                     'DD': {'d': 0.1273, 'k': 800}}

        # Convert k in kJ.mol-1.nm-2:
        for bond in bond_type.values():
            bond['k'] *= 4.184 * 100

        index_list = []
        for i in range(2, 8):
            index_list.append([1, i])
        bond_type['ID']['list'] = index_list

        index_list = [[2, 3], [2, 5], [2, 6], [2, 7],
                      [4, 3], [4, 5], [4, 6], [4, 7],
                      [5, 6], [5, 7], [3, 6], [3, 7]]
        # To avoid to much constraints
        # avoid with virtual sites ??
        # if vsite != 'none':
        #    index_list = []
        bond_type['DD']['list'] = index_list

        # Unit: angle in (degreee), k in kcal.mol-1.rad-2
        angle_type = {'DiMDi': {'theta': 180.0, 'k': 250},
                      'DiMDj': {'theta': 90.0, 'k': 250},
                      'MDiDj': {'theta': 45.0, 'k': 250},
                      'DiDjDi': {'theta': 90.0, 'k': 250},
                      'DiDjDk': {'theta': 60.0, 'k': 250}}

        # Convert k in kJ.mol-1.rad-2:
        for angle in angle_type.values():
            angle['k'] *= 4.184

        index_list = [[2, 1, 4], [3, 1, 5], [6, 1, 7]]
        angle_type['DiMDi']['list'] = index_list

        index_list = [[2, 1, 3], [2, 1, 5], [2, 1, 6],
                      [2, 1, 7], [4, 1, 3], [4, 1, 5],
                      [4, 1, 6], [4, 1, 7], [5, 1, 6],
                      [5, 1, 7], [3, 1, 6], [3, 1, 7]]
        angle_type['DiMDj']['list'] = index_list

        index_list = [[1, 2, 3], [1, 2, 5], [1, 2, 6],
                      [1, 2, 7], [1, 4, 3], [1, 4, 5],
                      [1, 4, 6], [1, 4, 7], [1, 5, 6],
                      [1, 5, 7], [1, 3, 6], [1, 3, 7]]
        angle_type['MDiDj']['list'] = index_list

        index_list = [[2, 3, 4], [2, 5, 4], [2, 6, 4],
                      [2, 7, 4], [3, 6, 5], [3, 7, 5],
                      [3, 2, 5], [3, 4, 5], [6, 2, 7],
                      [6, 3, 7], [6, 4, 7], [6, 5, 7]]
        angle_type['DiDjDi']['list'] = index_list

        index_list = [[2, 3, 6], [2, 3, 7], [2, 5, 6],
                      [2, 5, 7], [2, 6, 3], [2, 6, 5],
                      [2, 7, 3], [2, 7, 5], [4, 3, 6],
                      [4, 3, 7], [4, 5, 6], [4, 5, 7],
                      [4, 6, 3], [4, 6, 5], [4, 7, 3],
                      [4, 7, 5]]
        angle_type['DiDjDk']['list'] = index_list

        itp_file_name = self.top_file[:-4] + '_octa_ions.itp'
        # Create empty itp file:
        with open(itp_file_name, 'w') as filout:
            filout.write("; Itp file created by " + __name__ + "\n\n")

        ion_itp = Itp(name='Ion_octa_dummy',
                      fullname=itp_file_name.split('/')[-1],
                      path=itp_file_name)

        # Write the ion itp files:
        for name in ion_name:

            local_top = TopMol(name + 'D', 3)
            ion_itp.top_mol_list.append(local_top)

            # ATOM part:
            i = 1
            atom = {"num": i, "atom_type": name, "atom_name": name,
                    "res_num": 1, "res_name": name + 'D', "charge_num": i,
                    "charge": atomtypes[name]['charge'],
                    "mass": atomtypes[name]['mass']}
            local_top.atom_dict[i] = atom
            # Dummy atoms part:
            for i in range(2, 8):
                local_top.atom_dict[i] = {"num": i, "atom_type": 'D',
                                          "atom_name": 'D{}'.format(i - 1),
                                          "res_num": 1, "res_name": name + 'D',
                                          "charge_num": i,
                                          "charge": atomtypes['D']['charge'],
                                          "mass": atomtypes['D']['mass']}
            # BOND part:
            funct = 1
            for bond in bond_type.values():
                for i, j in bond['list']:
                    local_top.bond_list.append({'ai': i, 'aj': j,
                                                'funct': funct,
                                                'r': bond['d'],
                                                'k': bond['k']})
            # ANGLE part:
            funct = 1
            for angle in angle_type.values():
                for i, j, k in angle['list']:
                    local_top.angl_list.append({'ai': i, 'aj': j,
                                                'ak': k,
                                                'funct': funct,
                                                'theta': angle['theta'],
                                                'cth': angle['k']})

        ion_itp.write_file(itp_file_name)

        return(ion_itp)

    def switch_ion_octa_dummy(self, ion_name=['MN', 'ZN']):

        # Unit: A: kcal^(1/2).mol^(-1/2).Å^(-6) B: kcal^(1/2).mol^(-1/2).Å^(-3)
        atomtypes = {'NI': {'atnum': 28, 'mass': 40.69, 'charge': -1.0,
                            'A': 113.0, 'B': 84.0},
                     'CO': {'atnum': 27, 'mass': 40.93, 'charge': -1.0,
                            'A': 61.0, 'B': 31.0},
                     'ZN': {'atnum': 30, 'mass': 47.39, 'charge': -1.0,
                            'A': 68.0, 'B': 38.0},
                     'MN': {'atnum': 25, 'mass': 36.94, 'charge': -1.0,
                            'A': 171.0, 'B': 35.0},
                     'FE': {'atnum': 26, 'mass': 37.85, 'charge': -1.0,
                            'A': 70.0, 'B': 10.0},
                     'MG': {'atnum': 12, 'mass': 6.30, 'charge': -1.0,
                            'A': 63.0, 'B': 9.0},
                     'CA': {'atnum': 20, 'mass': 22.08, 'charge': -1.0,
                            'A': 350.0, 'B': 15.0},
                     'D': {'atnum': 0, 'mass': 3.0, 'charge': 0.5,
                           'A': 0.05, 'B': 0.0}}
        # Correct coodinates:
        ion_num = {}
        coor = pdb_manip.Coor(self.coor_file)
        for name in ion_name:
            print(name)
            ion_num[name] = coor.select_part_dict({'name': [name]}).num
            index = coor.get_index_selection({'res_name': [name]})
            print(name, len(index))
            coor.change_index_pdb_field(index,
                                        {'res_name': name + 'D'})
        coor.correct_ion_octa(ion_name)

        coor.write_pdb(self.coor_file[:-4] + '_newion.pdb',
                       check_file_out=False)
        self.coor_file = self.coor_file[:-4] + '_newion.pdb'

        # Correct topologie
        ion_atom_type = self.create_itp_atomtype_ion_octa_dummy(
            atomtypes, ion_name)
        ion_itp = self.create_itp_ion_octa_dummy(atomtypes, ion_name)

        # Create and add to topologie the atomtype itp:
        sys_top = TopSys(self.top_file)
        sys_top.add_atomtypes(ion_atom_type)

        # Use add_mol function !!!
        # But before need to use the remove mol function !!!
        sys_top.remove_ion(ion_name)

        for ion, num in ion_num.items():
            if num > 0:
                print(ion, sys_top.mol_num(ion))
                # sys_top.remove_mol(mol_name=ion)
                sys_top.add_mol(mol_name=ion+'D',
                                mol_itp_file=ion_itp.path,
                                mol_num=num)

        # sys_top.itp_list += [ion_itp]
        sys_top.write_file(self.top_file[:-4] + '_newion.top')
        self.top_file = self.top_file[:-4] + '_newion.top'

    @staticmethod
    def set_coor_aa_prot(coor_in, res_prot_dict, ff):
        """ Set manually residue protonation.

        :param coor_in: coordinate to update
        :type coor_in: Coor

        :param res_prot_dict: Dictionary of protonated residues
        :type res_prot_dict: dict

        :param ff: forcefield
        :type ff: str

        """

        res_dict = {'ASP':  ['ASP', 'ASPP', 'ASH'],
                    'ASPP': ['ASP', 'ASPP', 'ASH'],
                    'ASH':  ['ASP', 'ASPP', 'ASH'],
                    'GLUP': ['GLU', 'GLUP', 'GLH'],
                    'GLH':  ['GLU', 'GLUP', 'GLH'],
                    'GLU':  ['GLU', 'GLUP', 'GLH'],
                    'RN1': ['ARG', 'RN1'],
                    'LSN': ['LYS', 'LSN'],
                    'HSP': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'HSD': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'HSE': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'HIP': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'HID': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'HIE': ['HIS', 'HSP', 'HSD', 'HSE', 'HIP', 'HID', 'HIE'],
                    'CYM': ['CYS', 'CYM']}

        for res in res_prot_dict:

            local_dict = res_prot_dict[res]
            local_dict['res_name'] = res_dict[res]

            change_index_list = coor_in.get_index_selection(local_dict)

            if len(change_index_list) == 0:
                logger.warning("Selection dict {} return no"
                               " atoms. Check your selection".format(
                                   local_dict))
            else:
                logger.info("Set resname of sel dict {} to "
                            " {}.".format(
                                   local_dict, res))
                coor_in.change_index_pdb_field(index_list=change_index_list,
                                               change_dict={"res_name": res})

        chain_list = coor_in.get_attribute_selection(attribute='chain')
        res_prot_dict = {'ASP': [],
                         'ASPP': [],
                         'ASH': [],
                         'GLUP': [],
                         'GLH': [],
                         'GLU': [],
                         'RN1': [],
                         'LSN': [],
                         'HSP': [],
                         'HSD': [],
                         'HSE': [],
                         'HIP': [],
                         'HID': [],
                         'HIE': [],
                         'CYM': []}

        prot_result = {}

        logger.info('Protonation is :')
        for chain in chain_list:

            logger.info('Chain {}'.format(chain))

            prot_result[chain] = copy.deepcopy(res_prot_dict)

            chain_sel = coor_in.select_part_dict(selec_dict={'chain': [chain]})
            for resname in res_prot_dict.keys():
                local_sel = chain_sel.select_part_dict(
                    selec_dict={'res_name': [resname]})
                res_list = local_sel.get_attribute_selection(
                    attribute='res_num')
                if len(res_list) > 0:
                    prot_result[chain][resname] = res_list
                    logger.info('\tresidue {}: {}'.format(resname, res_list))

    def prepare_top_ligand(self, out_folder, name=None,
                           ff="amber99sb-ildn", water_model='tip3p',
                           include_mol={}):
        """Prepare the topologie of a ligand:

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ff: forcefield
        :type ff: str, optional, default="amber99sb-ildn"

        :param include_mol: list of ligand's residue name to include
        :type include_mol: list, optional, default=[]

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> # Create the topologie of a protein and do a minimisation:
        >>> lig = GmxSys(name='1D30', coor_file=TEST_PATH+'/1D30.pdb')



        .. note::
            Starting file need to be a pdb, this should be changed.
        """

        start_dir = os.path.abspath(".")

        # Create and go in out_folder:
        # This is necessary for the topologie creation
        os_command.create_and_go_dir(out_folder)

        # If name is not define use the object name
        if name is None:
            name = self.name

        # Save initial pdb file:
        start_pdb = self.coor_file
        start_coor = pdb_manip.Coor(start_pdb)

        # Get initial list of residue name
        res_name_input = start_coor.get_attribute_selection(
            attribute='res_name')

        mol_sys_list = []
        for resname in res_name_input:
            if resname in include_mol:
                mol_top = ambertools.make_amber_top_mol_rdkit(
                        start_pdb, resname, smile=include_mol[resname],
                        remove_h=False)
                mol_top['name'] = resname
                mol_sys_list.append(mol_top)
            else:
                logger.warning("residue(s) {} not included,"
                               " if you want add this residue, "
                               "add the residue name in include_mol".format(
                                    resname))

        # Create empty topologie:
        # Need to use write instead of append
        # In case the file aldready exists
        open('{}.top'.format(name), 'w').close()
        # Create empty coordinates:
        open('{}.pdb'.format(name), 'w').close()

        self.coor_file = '{}.pdb'.format(name)
        self.top_file = '{}.top'.format(name)

        # Get the system topologie:
        sys_top = TopSys('{}.top'.format(name))

        # Get path of forcefield and water model
        for forcefield in FORCEFIELD_PATH_LIST:
            if os_command.check_file_exist(
                        os.path.join(forcefield, ff + '.ff',
                                     'forcefield.itp')):
                path_ff = os.path.abspath(
                    os.path.join(forcefield, ff + '.ff',
                                 'forcefield.itp'))
            if os_command.check_file_exist(
                        os.path.join(forcefield, ff + '.ff',
                                     water_model + '.itp')):
                path_water = os.path.abspath(
                    os.path.join(forcefield, ff + '.ff',
                                 water_model + '.itp'))

        sys_top.forcefield = {'name': ff,
                              'fullname': "{}.ff/forcefield.itp".format(ff),
                              'path': path_ff}
        water_itp = Itp(name=water_model,
                        fullname="{}.ff/{}.itp".format(ff, water_model),
                        path=path_water)
        sys_top.itp_list = [water_itp]

        pdb_mol_list = []
        # Add the molecule in the sys topologie and update the water num:
        for mol in mol_sys_list:
            logger.info('Add Molecule: {}'.format(mol['name']))

            pdb_mol_list.append(mol['coor'])
            # Add topologie:
            # mol['GmxSys'].display()
            mol_top = TopSys(mol['GmxSys'].top_file)
            mol_itp = mol_top.get_include_no_posre_file_list()
            sys_top.add_mol(mol_name=mol['name'],
                            mol_itp_file=mol_itp[-1],
                            mol_num=mol['num'])
            # Add atomtypes itp:
            sys_top.add_atomtypes(mol_itp[0])

        if mol_sys_list:
            # Add coordinates:
            GmxSys.concat_coor(self.coor_file, *pdb_mol_list,
                               pdb_out=self.coor_file[:-4] + '_mol.pdb',
                               check_file_out=False)
            self.coor_file = self.coor_file[:-4] + '_mol.pdb'
            # Save topologie
            sys_top.write_file(self.top_file[:-4] + '_mol.top')
            self.top_file = self.top_file[:-4] + '_mol.top'

        os.chdir(start_dir)

    def cyclic_peptide_top(self, out_folder, name=None, check_file_out=True,
                           ff="charmm36-jul2017"):
        """Prepare a topologie for a cyclic peptide

            1. Create a peptide topologie with NH3+ Cter and COO- Nter using
            ``add_top()``.
            2. Delete useless termini atoms.
            3. Change atom types, names and charges.
            4. Add backbone bonds, angle and dihedral parameters.
            5. Finally compute the topologie with pdb2gmx add_top()

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> cyclic_pep = GmxSys(name='5vav', coor_file=TEST_PATH+'/5vav.pdb')
        >>>
        >>> #Basic usage :
        >>> cyclic_pep.cyclic_peptide_top(out_folder=os.path.join(str(\
TEST_OUT),'cyclic/top')) #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f ...test_files/5vav.pdb -o no_cyclic_5vav_pdb2gmx.pdb \
-p no_cyclic_5vav_pdb2gmx.top -i no_cyclic_5vav_posre.itp -water tip3p -ff \
charmm36-jul2017 -ignh -ter -vsite none
        Molecule topologie present in no_cyclic_5vav_pdb2gmx.top , \
extract the topologie in a separate file: no_cyclic_5vav_pdb2gmx.itp
        Protein_chain_A
        - ITP file: no_cyclic_5vav_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: no_cyclic_5vav_pdb2gmx.top
        Read rtp file : ...charmm36-jul2017.ff/merged.rtp
        Correct residue GLY  atom N    atom type NH3  to NH1...
        Correct residue GLY  atom HN   atom type HC   to H  ...
        Correct residue ASP  atom C    atom type CC   to C  ...
        Correct residue ASP  atom O    atom type OC   to O  ...
        Protein_chain_A
        Succeed to read file ...cyclic/top/no_cyclic_5vav_pdb2gmx.pdb ,  \
212 atoms found
        Succeed to save file ...cyclic/top/5vav_pdb2gmx.pdb
        >>> cyclic_top = TopSys(cyclic_pep.top_file)
        >>> print(cyclic_top.charge())
        0.0
        >>> cyclic_top.prot_res_num()
        Protein_chain_A : 1
        Get Res num of Protein_chain_A : 14
        Total number of residue: 14
        14
        >>> cyclic_top.display() #doctest: +ELLIPSIS
        Forcefield include :
         charmm36-jul2017
        - ITP file: 5vav_pdb2gmx
        - molecules defined in the itp file:
        * Protein_chain_A
        - ITP file: tip3p
        - molecules defined in the itp file:
        * SOL
        - ITP file: ions
        - molecules defined in the itp file:
        * OH
        * LI
        * NA
        * K
        * CS
        * CL
        * CA
        * MG
        * ZN
        Mol List:
           * 1 Protein_chain_A
        Mol Name:
         CYC-MC12
        >>> cyclic_pep.em(out_folder=TEST_OUT+'/cyclic/em/', nsteps=10, \
create_box_flag=True) #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../cyclic/top/5vav_pdb2gmx.pdb -o \
.../cyclic/top/5vav_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        - Create the tpr file 5vav.tpr
        gmx grompp -f 5vav.mdp -c ../top/5vav_pdb2gmx_box.pdb -r \
../top/5vav_pdb2gmx_box.pdb -p ../top/5vav_pdb2gmx.top -po out_5vav.mdp \
-o 5vav.tpr -maxwarn 1
        - Launch the simulation 5vav.tpr
        gmx mdrun -s 5vav.tpr -deffnm 5vav -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        >>> cyclic_amber_pep = GmxSys(name='5vav_amber', \
coor_file=TEST_PATH+'/5vav.pdb')
        >>> cyclic_amber_pep.cyclic_peptide_top(out_folder=\
os.path.join(str(TEST_OUT),'cyclic/top'),ff='amber99sb-ildn')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/5vav.pdb -o \
no_cyclic_5vav_amber_pdb2gmx.pdb -p no_cyclic_5vav_amber_pdb2gmx.top \
-i no_cyclic_5vav_amber_posre.itp -water tip3p -ff amber99sb-ildn -ignh \
-ter -vsite none
        Molecule topologie present in no_cyclic_5vav_amber_pdb2gmx.top , \
extract the topologie in a separate file: no_cyclic_5vav_amber_pdb2gmx.itp
        Protein_chain_A
        - ITP file: no_cyclic_5vav_amber_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: no_cyclic_5vav_amber_pdb2gmx.top
        Read rtp file : ...amber99sb-ildn.ff/aminoacids.rtp
        Correct residue GLY  atom N    atom type N3   to N ...
        Correct residue GLY  atom HA1  atom type HP   to H1...
        Correct residue GLY  atom HA2  atom type HP   to H1...
        Correct residue ASP  atom O    atom type O2   to O ...
        Protein_chain_A
        Succeed to read file ...cyclic/top/no_cyclic_5vav_amber_pdb2gmx.pdb \
,  212 atoms found
        Succeed to save file ...cyclic/top/5vav_amber_pdb2gmx.pdb
        >>> cyclic_amber_top = TopSys(cyclic_amber_pep.top_file)
        >>> print(cyclic_amber_top.charge())
        0.0
        >>> cyclic_amber_pep.em(out_folder=TEST_OUT+'/cyclic/em/', \
nsteps=10, create_box_flag=True) #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f ...cyclic/top/5vav_amber_pdb2gmx.pdb -o \
...cyclic/top/5vav_amber_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        - Create the tpr file 5vav_amber.tpr
        gmx grompp -f 5vav_amber.mdp -c ../top/5vav_amber_pdb2gmx_box.pdb \
-r ../top/5vav_amber_pdb2gmx_box.pdb -p ../top/5vav_amber_pdb2gmx.top -po \
out_5vav_amber.mdp -o 5vav_amber.tpr -maxwarn 1
        - Launch the simulation 5vav_amber.tpr
        gmx mdrun -s 5vav_amber.tpr -deffnm 5vav_amber -nt 0 -ntmpi 0 \
-nsteps -2 -nocopyright

        .. note::
            No options are allowed (water model, termini capping) except
            for vsites.
        .. warning::
            Has not been tested with special residues like GLY or PRO !!

        """

        N_ter_dic = {"NH3+": "0", "NH2": "1", "5TER": "2", "None": "3"}
        C_ter_dic = {"COO-": "0", "COOH": "1", "CT2": "2", "3TER": "3",
                     "None": "4"}

        # If name is not define use the object name
        if name is None:
            name = self.name

        if check_file_out and os_command.check_file_and_create_path(
                os.path.join(out_folder, name + "_pdb2gmx.top")):
            logger.info("create_top not launched {} already exist".format(
                out_folder + "/" + name + "_pdb2gmx.top"))
            self.coor_file = os.path.join(out_folder, name + "_pdb2gmx.pdb")
            self.top_file = os.path.join(out_folder, name + "_pdb2gmx.top")
            return

        # Create peptide topologie with NH2 Cter and COO- Nter
        self.add_top(out_folder=out_folder, name="no_cyclic_" + name,
                     water="tip3p", ff=ff,
                     pdb2gmx_option_dict={'vsite': 'none',
                                          'ignh': None,
                                          'ter': None},
                     check_file_out=False,
                     input_pdb2gmx=N_ter_dic["NH3+"] + "\n" +
                     C_ter_dic["COO-"])

        # Make the top clean:
        top_pep = TopSys(self.top_file)
        mol_top = top_pep.itp_list[0].top_mol_list[0]
        res_num = mol_top.get_res_num()

        # Define termini atoms:
        # C-ter NH3

        if ff.startswith('amber'):
            # N_type = 'N'
            # HN_type = 'H'
            HN_name = 'H'
            # GLY_HA_type = 'H1'
            angle_func = 1
            dihe_func = 9
            dihe_impr_func = 4
        elif ff.startswith('charmm'):
            # For charmm
            # N_type = 'NH1'
            # HN_type = 'H'
            HN_name = 'HN'
            angle_func = 5
            dihe_func = 9
            dihe_impr_func = 2

        # Delete useless ter atoms (HN1 and HN2 are for PRO):
        to_del_name_n_ter = ['H2', 'H3', 'HN2', 'HN1']
        del_index = mol_top.get_selection_index(
            selec_dict={'atom_name': to_del_name_n_ter, 'res_num': [1]}) +\
            mol_top.get_selection_index(
                selec_dict={'atom_name': ['OT2', 'OC2'], 'res_num': [res_num]})

        mol_top.delete_atom(index_list=del_index)

        # Change atom name:
        #
        # chg_index = mol_top.get_selection_index(
        #   selec_dict={'atom_name': ['N'], 'res_num': [1]})[0]
        # mol_top.atom_dict[chg_index]['atom_type'] = N_type
        # mol_top.atom_dict[chg_index]['charge'] = -0.470

        if mol_top.atom_dict[1]['res_name'] != 'PRO':
            chg_index = mol_top.get_selection_index(
                selec_dict={'atom_name': ['H1', 'HN1'], 'res_num': [1]})[0]
            mol_top.atom_dict[chg_index]['atom_name'] = HN_name

        # mol_top.atom_dict[chg_index]['atom_type'] = HN_type
        # mol_top.atom_dict[chg_index]['charge'] = 0.310
        # chg_index = mol_top.get_selection_index(
        #   selec_dict={'atom_name': ['CA'], 'res_num': [1]})[0]
        # mol_top.atom_dict[chg_index]['charge'] = \
        #   mol_top.atom_dict[chg_index]['charge'] - 0.12
        chg_index = mol_top.get_selection_index(
            selec_dict={'atom_name': ['OT1', 'OC1'], 'res_num': [res_num]})[0]
        # mol_top.atom_dict[chg_index]['atom_type'] = 'O'
        mol_top.atom_dict[chg_index]['atom_name'] = 'O'
        # mol_top.atom_dict[chg_index]['charge'] = -0.51

        # if (mol_top.atom_dict[1]['res_name'] == 'GLY' and
        #       ff.startswith('amber')):
        #    chg_index = mol_top.get_selection_index(
        #       selec_dict={'atom_name': ['HA1'], 'res_num': [1]})[0]
        #    mol_top.atom_dict[chg_index]['atom_type'] = GLY_HA_type
        #    chg_index = mol_top.get_selection_index(
        #       selec_dict={'atom_name': ['HA2'], 'res_num': [1]})[0]
        #    mol_top.atom_dict[chg_index]['atom_type'] = GLY_HA_type

        # chg_index = mol_top.get_selection_index(
        #   selec_dict={'atom_name': ['C'], 'res_num': [res_num]})[0]
        # mol_top.atom_dict[chg_index]['atom_type'] = 'C'
        # mol_top.atom_dict[chg_index]['charge'] = 0.51

        last_res_index = chg_index

        # Get index for residue i-2
        prev_2_C_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['C'], 'res_num': [res_num - 1]})[0]
        # Get index for residue i-1
        prev_C_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['C'], 'res_num': [res_num]})[0]
        prev_O_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['O'], 'res_num': [res_num]})[0]
        prev_CA_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['CA'], 'res_num': [res_num]})[0]
        prev_N_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['N'], 'res_num': [res_num]})[0]

        # check if res is GLY:
        if mol_top.atom_dict[last_res_index]['res_name'] != 'GLY':
            prev_HA_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA'], 'res_num': [res_num]})[0]
            prev_CB_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['CB'], 'res_num': [res_num]})[0]
        else:
            prev_HA_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA1'], 'res_num': [res_num]})[0]
            prev_CB_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA2'], 'res_num': [res_num]})[0]

        # Get index for residue i
        C_index = mol_top.get_selection_index(
            selec_dict={'atom_name': ['C'], 'res_num': [1]})[0]
        # O_index = mol_top.get_selection_index(
        #   selec_dict={'atom_name': ['O'], 'res_num': [1]})[0]
        CA_index = mol_top.get_selection_index(
            selec_dict={'atom_name': ['CA'], 'res_num': [1]})[0]
        N_index = mol_top.get_selection_index(
            selec_dict={'atom_name': ['N'], 'res_num': [1]})[0]
        # check if res is PRO:
        if mol_top.atom_dict[1]['res_name'] != 'PRO':
            HN_index = mol_top.get_selection_index(selec_dict={
                'atom_name': [HN_name], 'res_num': [1]})[0]
        else:
            HN_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['CD'], 'res_num': [1]})[0]

        # check if res is GLY:
        if mol_top.atom_dict[1]['res_name'] != 'GLY':
            HA_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA'], 'res_num': [1]})[0]
            CB_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['CB'], 'res_num': [1]})[0]
        else:
            HA_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA1'], 'res_num': [1]})[0]
            CB_index = mol_top.get_selection_index(selec_dict={
                'atom_name': ['HA2'], 'res_num': [1]})[0]
        # Get index for residue i+1
        next_N_index = mol_top.get_selection_index(selec_dict={
            'atom_name': ['N'], 'res_num': [2]})[0]

        # Add backbone bonds, angle, dihedral parameters:
        # Bond:
        # N-C
        mol_top.bond_list.append({'ai': N_index, 'aj': prev_C_index,
                                  'funct': 1, 'r': '', 'k': ''})

        # Pairs
        pair_list = [[N_index, prev_HA_index],
                     [N_index, prev_N_index],
                     [N_index, prev_CB_index],
                     [HN_index, prev_O_index],
                     [HN_index, prev_CA_index],
                     [CA_index, prev_O_index],
                     [CA_index, prev_CA_index],
                     [HA_index, prev_C_index],
                     [CB_index, prev_C_index]]

        for ai, aj in pair_list:
            mol_top.pair_list.append({'ai': ai, 'aj': aj, 'funct': 1})

        # Angle:

        angle_list = [[N_index, prev_C_index, prev_O_index],
                      [N_index, prev_C_index, prev_CA_index],
                      [HN_index, N_index, prev_C_index],
                      [CA_index, N_index, prev_C_index]]

        for ai, aj, ak in angle_list:
            mol_top.angl_list.append({'ai': ai, 'aj': aj, 'ak': ak,
                                      'funct': angle_func,
                                      'theta': '', 'cth': ''})

        # Dihed: type 9
        dihed_list = [[N_index, prev_C_index, prev_CA_index, prev_HA_index],
                      [N_index, prev_C_index, prev_CA_index, prev_CB_index],
                      [N_index, prev_C_index, prev_CA_index, prev_N_index],
                      [HN_index, N_index, prev_C_index, prev_O_index],
                      [HN_index, N_index, prev_C_index, prev_CA_index],
                      [CA_index, N_index, prev_C_index, prev_O_index],
                      [CA_index, N_index, prev_C_index, prev_CA_index],
                      [HA_index, CA_index, N_index, prev_C_index],
                      [CB_index, CA_index, N_index, prev_C_index],
                      [C_index, CA_index, N_index, prev_C_index]]

        for ai, aj, ak, al in dihed_list:
            mol_top.dihe_list.append({'ai': ai, 'aj': aj, 'ak': ak,
                                      'al': al, 'funct': dihe_func,
                                      'phase': '', 'kd': '', 'pn': ''})

        # Dihed: type 2
        if ff.startswith('charmm'):
            dihed_list_impr = [[prev_C_index, prev_CA_index, N_index,
                                prev_O_index],
                               [N_index, prev_C_index, CA_index, HN_index]]
        elif ff.startswith('amber'):
            dihed_list_impr = [[prev_CA_index, N_index, prev_C_index,
                                prev_O_index],
                               [prev_C_index, CA_index, N_index, HN_index]]

        for ai, aj, ak, al in dihed_list_impr:
            mol_top.dihe_list.append({'ai': ai, 'aj': aj, 'ak': ak,
                                      'al': al, 'funct': dihe_impr_func,
                                      'phase': '', 'kd': '', 'pn': ''})

        # Cmap
        if ff.startswith('charmm'):
            cmap_list = [[prev_2_C_index, prev_N_index, prev_CA_index,
                          prev_C_index, N_index],
                         [prev_C_index, N_index, CA_index, C_index,
                          next_N_index]]

            for ai, aj, ak, al, am in cmap_list:
                mol_top.cmap_list.append(
                    {'ai': ai, 'aj': aj, 'ak': ak,
                     'al': al, 'am': am, 'funct': 1})

        # Correct charge and atom type base on ff .rtp file
        mol_top.correct_charge_type(forcefield=top_pep.forcefield)

        # Save itp:

        top_pep.itp_list[0].write_file(os.path.join(
            out_folder, name + "_pdb2gmx.itp"))
        top_pep.itp_list[0].name = name + "_pdb2gmx.itp"
        top_pep.itp_list[0].fullname = name + "_pdb2gmx.itp"
        top_pep.itp_list[0].path = os.path.abspath(
            os.path.join(out_folder, name + "_pdb2gmx.itp"))
        # Save top:
        top_pep.write_file(os.path.join(out_folder, name + "_pdb2gmx.top"))
        self.top_file = os.path.join(out_folder, name + "_pdb2gmx.top")

        # Correct pdb file:
        coor_pep = pdb_manip.Coor(self.coor_file)
        to_del_index = coor_pep.get_index_selection(
            selec_dict={'name': to_del_name_n_ter, 'res_num': [1]})\
            + coor_pep.get_index_selection(
                selec_dict={'name': ['OT2', 'OC2'], 'res_num': [res_num]})

        coor_pep.del_atom_index(to_del_index)

        # Change atom name :

        chg_index = coor_pep.get_index_selection(
            selec_dict={'name': ['H1'], 'res_num': [1]})
        coor_pep.change_index_pdb_field(chg_index, {'name': HN_name})

        chg_index = coor_pep.get_index_selection(
            selec_dict={'name': ['OT1', 'OC1'], 'res_num': [res_num]})
        coor_pep.change_index_pdb_field(chg_index, {'name': 'O'})

        # save pdb:
        coor_pep.write_pdb(pdb_out=os.path.join(out_folder,
                                                name + "_pdb2gmx.pdb"))
        self.coor_file = os.path.join(out_folder, name + "_pdb2gmx.pdb")

        # First read and save to fix the molecule top include in .top:
        # top = TopSys(self.top_file)
        # top.write_file(self.top_file)

        # Fix the molecule posre files with the wrong atom number in the .top:

        top_pep.add_posre(posre_name="POSRES", selec_dict={
            'atom_name': HA_NAME},
            fc=[1000, 1000, 1000])
        top_pep.add_posre(posre_name="POSRES_HA_LOW", selec_dict={
            'atom_name': HA_NAME},
            fc=[100, 100, 100])
        top_pep.add_posre(posre_name="POSRES_CA", selec_dict={
            'atom_name': ['CA']},
            fc=[1000, 1000, 1000])
        top_pep.add_posre(posre_name="POSRES_CA_LOW", selec_dict={
            'atom_name': ['CA']},
            fc=[100, 100, 100])

    def add_disulfide_bonds(self, res_list, out_folder, name=None,
                            check_file_out=True, ff="charmm36-jul2017"):
        """Add disulfide bonds to a single protein topologie.
        Topologie has to be computed before using this function.
        Set specifically which cystein residues need to be bonded:
        Example of res_list = [[4, 7], [10, 20]], will connect
        cystein residues 4 to 7 and 10 to 20.

        :param res_list: list of list of cystein residues to be bonded
        :type res_listr: list

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param ff: forcefield.
        :type ff: str, optional, default="charmm36-jul2017"

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        Add the following terms:

        * 1 Bond

            * SG-SG

        * 2 Angle

            * SG-SG-CB
            * CB-SG-SG

        * 7 Dihed

            * CA-CB-SG-SG
            * HB1-CB-SG-SG
            * HB2-CB-SG-SG
            * CB-SG-SG-CB
            * SG-SG-CB-HB1
            * SG-SG-CB-HB2
            * SG-SG-CB-CA


        :Example:

        >>> show_log()
        >>> TEST_OUT = getfixture('tmpdir')
        >>>
        >>> # Measure s-s bond length:
        >>> ss_coor = pdb_manip.Coor(TEST_PATH+'/1dn3_cys.pdb')  \
#doctest: +ELLIPSIS
        Succeed to read file .../1dn3_cys.pdb ,  \
144 atoms found
        >>> cystein_s_index = ss_coor.get_index_selection({'name': ['SG'], \
'res_name' : ['CYS']})
        >>> print(cystein_s_index)
        [85, 118]
        >>> distance = pdb_manip.Coor.atom_dist(ss_coor.atom_dict[\
cystein_s_index[0]], ss_coor.atom_dict[cystein_s_index[1]])
        >>> print('S-S distance = {:.2f} Å'.format(distance))
        S-S distance = 6.38 Å
        >>> no_ss_pep = GmxSys(name='1dn3_cys', coor_file=TEST_PATH+\
'/1dn3_cys.pdb')
        >>>
        >>> #Basic usage :
        >>> no_ss_pep.prepare_top(out_folder=os.path.join(str(\
TEST_OUT),'1dn3/top')) #doctest: +ELLIPSIS
        Succeed to read file .../1dn3_cys.pdb ,  144 atoms found
        Succeed to read file .../1dn3_cys.pdb ,  144 atoms found
        Succeed to save file tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method\
=propka --with-ph=7.00 tmp_pdb2pqr.pdb 00_1dn3_cys.pqr
        Succeed to read file 00_1dn3_cys.pqr ,  231 atoms found
        Chain: A  Residue: 0 to 14
        Succeed to save file 01_1dn3_cys_good_his.pdb
        - Create topologie
        gmx pdb2gmx -f 01_1dn3_cys_good_his.pdb -o 1dn3_cys_pdb2gmx.pdb \
-p 1dn3_cys_pdb2gmx.top -i 1dn3_cys_posre.itp -water tip3p -ff \
charmm36-jul2017 -ignh -vsite none
        Molecule topologie present in 1dn3_cys_pdb2gmx.top , extract the \
topologie in a separate file: 1dn3_cys_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1dn3_cys_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1dn3_cys_pdb2gmx.top
        >>> no_ss_pep.add_disulfide_bonds(res_list=[[9, 12]], \
out_folder=os.path.join(str(TEST_OUT),'1dn3/top_ss')) #doctest: +ELLIPSIS
        Succeed to read file .../1dn3/top/1dn3_cys_pdb2gmx.pdb ,  \
231 atoms found
        Succeed to save file .../1dn3/top_ss/1dn3_cys_ss_bond.pdb
        Read rtp file : .../charmm36-jul2017.ff/merged.rtp
        Correct residue CYS2 atom SG   atom type S    to SM ...
        Correct residue CYS2 atom SG   atom type S    to SM ...
        Protein_chain_A
        >>> no_ss_pep.em(out_folder=TEST_OUT+'/1dn3/em_ss/', \
nsteps=10, create_box_flag=True) #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f ...1dn3/top_ss/1dn3_cys_ss_bond.pdb -o \
.../1dn3/top_ss/1dn3_cys_ss_bond_box.pdb -bt dodecahedron -d 1.0
        - Create the tpr file 1dn3_cys.tpr
        gmx grompp -f 1dn3_cys.mdp -c ../top_ss/1dn3_cys_ss_bond_box.pdb -r \
../top_ss/1dn3_cys_ss_bond_box.pdb -p ../top_ss/1dn3_cys_ss_bond.top -po \
out_1dn3_cys.mdp -o 1dn3_cys.tpr -maxwarn 1
        - Launch the simulation 1dn3_cys.tpr
        gmx mdrun -s 1dn3_cys.tpr -deffnm 1dn3_cys -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        >>> # Need to convert the gro to pdb:
        >>> no_ss_pep.convert_trj(traj=False) #doctest: +ELLIPSIS
        - Convert trj/coor
        gmx trjconv -f .../em_ss/1dn3_cys.gro -o \
.../em_ss/1dn3_cys_compact.pdb -s .../em_ss/1dn3_cys.tpr -ur compact -pbc mol
        >>> # Measure s-s bond length:
        >>> ss_coor = pdb_manip.Coor(no_ss_pep.coor_file) #doctest: +ELLIPSIS
        Succeed to read file .../em_ss/1dn3_cys_compact.pdb ,  229 atoms found
        >>> cystein_s_index = ss_coor.get_index_selection({'name': ['SG'], \
'res_name' : ['CYS']})
        >>> print(cystein_s_index)
        [135, 189]
        >>> distance = pdb_manip.Coor.atom_dist(ss_coor.atom_dict[\
cystein_s_index[0]], ss_coor.atom_dict[cystein_s_index[1]])
        >>> print('S-S distance = {:.2f} Å'.format(distance)) \
#doctest: +ELLIPSIS
        S-S distance = 2.0... Å

        .. note::
            No options are allowed (water model, termini capping) except
            for vsites.

        """

        # If name is not define use the object name
        if name is None:
            name = self.name

        if check_file_out and os_command.check_file_and_create_path(
                os.path.join(out_folder, name + "_ss_bond.top")):
            logger.info("create_top not launched {} already exist".format(
                out_folder + "/" + name + "_pdb2gmx.top"))
            self.top_file = os.path.join(out_folder, name + "_ss_bond.top")
            return

        top_prot = TopSys(self.top_file)
        mol_top = top_prot.itp_list[0].top_mol_list[0]

        # Get cystein index:
        coor_prot = pdb_manip.Coor(self.coor_file)

        cys_S_index_list = []
        cys_H_index_list = []

        for res_pair in res_list:
            cys_pair_index = mol_top.get_selection_index({'atom_name': ['SG'],
                                                          'res_name': ['CYS'],
                                                          'res_num': res_pair})
            if len(cys_pair_index) != 2:
                logger.warning("Wrong residue number {}, at least"
                               " one residue is not a cystein".format(
                                res_pair))

            cys_S_index_list.append(cys_pair_index)
            cys_H_index_list += mol_top.get_selection_index({
                'atom_name': ['HG1'], 'res_name': ['CYS'],
                'res_num': res_pair})
        # print('S atoms', cys_S_index_list)
        # print('SH atoms', cys_H_index_list)

        if ff.startswith('amber'):
            bond_func = 1
            angle_func = 1
            dihe_func = 9
        elif ff.startswith('charmm'):
            bond_func = 1
            angle_func = 5
            dihe_func = 9

        # Delete SH atoms
        mol_top.delete_atom(index_list=cys_H_index_list)
        flat_res_list = [item for sublist in res_list for item in sublist]
        # print(flat_res_list)
        coor_cys_H_index_list = coor_prot.get_index_selection({
            'name': ['HG1'], 'res_name': ['CYS'], 'res_num': flat_res_list})
        coor_prot.del_atom_index(index_list=coor_cys_H_index_list)

        # Save coor
        coor_prot.write_pdb(os.path.join(out_folder, name + "_ss_bond.pdb"))
        self.coor_file = os.path.join(out_folder, name + "_ss_bond.pdb")

        for res_pair in res_list:

            SG_1_index = mol_top.get_selection_index({
                'atom_name': ['SG'], 'res_name': ['CYS'],
                'res_num': [res_pair[0]]})[0]
            SG_2_index = mol_top.get_selection_index({
                'atom_name': ['SG'], 'res_name': ['CYS'],
                'res_num': [res_pair[1]]})[0]
            # print('SG1', mol_top.atom_dict[SG_1_index])
            # print('SG2', mol_top.atom_dict[SG_2_index])

            CA_1_index = mol_top.get_selection_index({
                'atom_name': ['CA'], 'res_name': ['CYS'],
                'res_num': [res_pair[0]]})[0]
            CA_2_index = mol_top.get_selection_index({
                'atom_name': ['CA'], 'res_name': ['CYS'],
                'res_num': [res_pair[1]]})[0]
            CB_1_index = mol_top.get_selection_index({
                'atom_name': ['CB'], 'res_name': ['CYS'],
                'res_num': [res_pair[0]]})[0]
            CB_2_index = mol_top.get_selection_index({
                'atom_name': ['CB'], 'res_name': ['CYS'],
                'res_num': [res_pair[1]]})[0]
            HB1_1_index = mol_top.get_selection_index({
                'atom_name': ['HB1'], 'res_name': ['CYS'],
                'res_num': [res_pair[0]]})[0]
            HB1_2_index = mol_top.get_selection_index({
                'atom_name': ['HB1'], 'res_name': ['CYS'],
                'res_num': [res_pair[1]]})[0]
            HB2_1_index = mol_top.get_selection_index({
                'atom_name': ['HB2'], 'res_name': ['CYS'],
                'res_num': [res_pair[0]]})[0]
            HB2_2_index = mol_top.get_selection_index({
                'atom_name': ['HB2'], 'res_name': ['CYS'],
                'res_num': [res_pair[1]]})[0]

            # Add S-S bond
            mol_top.bond_list.append({'ai': SG_1_index,
                                      'aj': SG_2_index,
                                      'funct': bond_func,
                                      'r': '',
                                      'k': ''})

            # Add X-S-S, S-S-X angles
            local_angle_list = [[SG_1_index, SG_2_index, CB_2_index]]
            local_angle_list.append([CB_1_index, SG_1_index, SG_2_index])

            for angle in local_angle_list:
                mol_top.angl_list.append({'ai': angle[0],
                                          'aj': angle[1],
                                          'ak': angle[2],
                                          'funct': angle_func,
                                          'theta': '',
                                          'cth': ''})

            # Add X-S-S-X, S-S-X-X, X-X-S-S angles
            local_dihed_list = [[CA_1_index, CB_1_index, SG_1_index,
                                 SG_2_index]]
            local_dihed_list.append([HB1_1_index, CB_1_index, SG_1_index,
                                     SG_2_index])
            local_dihed_list.append([HB2_1_index, CB_1_index, SG_1_index,
                                     SG_2_index])
            local_dihed_list.append([CB_1_index, SG_1_index, SG_2_index,
                                     CB_2_index])
            local_dihed_list.append([SG_1_index, SG_2_index, CB_2_index,
                                     CA_2_index])
            local_dihed_list.append([SG_1_index, SG_2_index, CB_2_index,
                                     HB1_2_index])
            local_dihed_list.append([SG_1_index, SG_2_index, CB_2_index,
                                     HB2_2_index])

            for dihed in local_dihed_list:
                mol_top.dihe_list.append({'ai': dihed[0],
                                          'aj': dihed[1],
                                          'ak': dihed[2],
                                          'al': dihed[3],
                                          'funct': dihe_func,
                                          'phase': '',
                                          'kd': '',
                                          'pn': ''})

        # print(top_prot.itp_list[0].top_mol_list[0].atom_dict)

        # Correct charge and atom type base on ff .rtp file
        cys_index_list = mol_top.get_selection_index({
            'res_name': ['CYS'], 'res_num': flat_res_list})
        mol_top.correct_charge_type(forcefield=top_prot.forcefield,
                                    index_list=cys_index_list)

        # Save itp:
        # print('top name', top_prot.itp_list[0].name)
        top_prot.itp_list[0].write_file(os.path.join(
            out_folder, name + "_ss_bond.itp"))
        top_prot.itp_list[0].name = name + "_ss_bond.itp"
        top_prot.itp_list[0].fullname = name + "_ss_bond.itp"
        top_prot.itp_list[0].path = os.path.abspath(
            os.path.join(out_folder, name + "_ss_bond.itp"))

        # Save top:
        top_prot.write_file(os.path.join(out_folder, name + "_ss_bond.top"))
        self.top_file = os.path.join(out_folder, name + "_ss_bond.top")

        # Fix the molecule posre files with the wrong atom number in the .top:

        top_prot.add_posre(posre_name="POSRES", selec_dict={
            'atom_name': HA_NAME},
            fc=[1000, 1000, 1000])
        top_prot.add_posre(posre_name="POSRES_HA_LOW", selec_dict={
            'atom_name': HA_NAME},
            fc=[100, 100, 100])
        top_prot.add_posre(posre_name="POSRES_CA_LOW", selec_dict={
            'atom_name': ['CA']},
            fc=[100, 100, 100])
        top_prot.add_posre(posre_name="POSRES_CA", selec_dict={
            'atom_name': ['CA']},
            fc=[1000, 1000, 1000])

    #######################################################
    # ###########  SYSTEM CREATION FUNCTIONS  #############
    #######################################################

    def create_box(self, name=None, dist=1.0, box_type="dodecahedron",
                   check_file_out=True):
        """Create pbc box using ``gmx editconf``

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param dist: Distance between the solute and the box (nm)
        :type dist: float, default=1.0

        :param box_type: Box type ("triclinic", "cubic", "dodecahedron",
           "octahedron")
        :type box_type: str, default="dodecahedron"

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/create_box/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../create_box/top_SH3/1y0m_pdb2gmx.pdb -o \
.../create_box/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0

        .. note::
            If ``name`` is not defined, the command will create a new pdb file
            name after the input one and adding "_box.pdb".
            If ``name`` is defined the pdb filed will be saved in the same
            directory as input file, the "_box.pdb" will be added to ``name``.
        """

        logger.info("- Create pbc box")

        # If name is not define use the object coor name and add _box.pdb
        if name is None:
            box_coor = self.coor_file[:-4] + "_box.pdb"
        # If not use the coor_path Path and add name and _box.pdb
        else:
            box_coor = os.path.join(
                os_command.get_directory(self.coor_file), name + "_box.pdb")

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(box_coor):
            logger.info("create_box not launched {} already exist".format(
                box_coor))
            self.coor_file = box_coor
            return

        cmd_box = os_command.Command([GMX_BIN, "editconf",
                                      "-f", self.coor_file,
                                      "-o", box_coor,
                                      "-bt", box_type,
                                      "-d", str(dist)])

        cmd_box.display()
        cmd_box.run()

        self.coor_file = box_coor

    def convert_trj(self, name=None, ur="compact", pbc="mol", select="System",
                    traj=True, specific_coor_out=None, check_file_out=True,
                    tpr=None, **cmd_args):
        """Convert a trajectory or coordinate file using the commande
        ``gmx trjconv``.

        This is specially usefull when the protein is break across pbc. Using
        ``convert_trj()`` with default parameters will fix it.

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ur: unit-cell representation ("rect", "tric", "compact")
        :type ur: str, default="compact"

        :param pbc: PBC treatment ("none", "mol", "res",
           "atom", "nojump", "cluster", "whole")
        :type pbc: str, default="mol"

        :param select: group for output
        :type select: str, default="System"

        :param specific_coor_out: specific output file
        :type specific_coor_out: str, optional, default=None

        :param traj: Flag to convert trajectory or coordinates
        :type traj: bool, default=True

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param cmd_args: Optional arguments for ``gmx trjconv``

        **Object requirement(s):**

            * self.tpr
            * self.coor_file or self.xtc

        **Object field(s) changed:**

            * self.coor_file or self.xtc

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/convert_trj/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../convert_trj/top_SH3/1y0m_pdb2gmx.pdb -o \
.../convert_trj/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        >>> prot.solvate_box(out_folder=TEST_OUT+'/convert_trj/top_SH3_water/')
        - Solvate the pbc box
        Copy topologie file and dependancies
        >>> prot.em(out_folder=TEST_OUT+'/convert_trj/em_SH3_water/', \
nsteps=10, constraints="none")
        - Create the tpr file 1y0m.tpr
        gmx grompp -f 1y0m.mdp -c ../top_SH3_water/1y0m_water.pdb -r \
../top_SH3_water/1y0m_water.pdb -p ../top_SH3_water/1y0m_water.top -po \
out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
        - Launch the simulation 1y0m.tpr
        gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright
        >>> prot.convert_trj(traj=False) #doctest: +ELLIPSIS
        - Convert trj/coor
        gmx trjconv -f .../convert_trj/em_SH3_water/1y0m.gro -o \
.../convert_trj/em_SH3_water/1y0m_compact.pdb -s \
.../convert_trj/em_SH3_water/1y0m.tpr -ur compact -pbc mol

        .. note::
            If ``name`` is not defined, the command will create a new pdb
            file name after the input one and adding "_compact.pdb" or
            "_compact.xtc". If ``name`` is defined the pdb filed will be
            saved in the same directory as input file, the "_compact.pdb"
            or "_compact.xtc" will be added to ``name``.
        """

        logger.info("- Convert trj/coor")
        if traj:
            coor_in = self.xtc
            if name is None:
                coor_out = self.xtc[:-4] + "_compact.xtc"
            else:
                coor_out = os.path.join(os_command.get_directory(
                    self.xtc), name + "_compact.xtc")
        else:
            coor_in = self.coor_file
            if name is None:
                coor_out = self.coor_file[:-4] + "_compact.pdb"
            else:
                coor_out = os.path.join(os_command.get_directory(
                    self.coor_file), name + "_compact.pdb")

        if tpr is None:
            tpr = self.tpr

        if specific_coor_out is not None:
            coor_out = specific_coor_out

        # Check if output files exist:
        if check_file_out and os.path.isfile(coor_out):
            logger.info("convert trj not launched {} already exist".format(
                coor_out))
            if traj:
                self.xtc = coor_out
            else:
                self.coor_file = coor_out
            return

        if self.tpr is None:
            logger.info("tpr file missing, function \"convert_trj\" could "
                        "not be executed")
            raise ValueError("tpr file is missing")

        if self.ndx is not None:
            cmd_args.update({'n': self.ndx})

        cmd_convert = os_command.Command([GMX_BIN, "trjconv",
                                          "-f", coor_in,
                                          "-o", coor_out,
                                          "-s", tpr,
                                          "-ur", ur,
                                          "-pbc", pbc],
                                         **cmd_args)

        cmd_convert.display()
        cmd_convert.run(com_input=select)

        if traj:
            self.xtc = coor_out
        else:
            self.coor_file = coor_out

    def copy_box(self, nbox, name=None, check_file_out=True,
                 renumber=False, **cmd_args):
        """Copy images of a given corrdinates in x, y, and z directions using
        ``gmx genconf``.

        nbox needs a list of 3 string for number x,y,z dimensions copy

        This is specially usefull when the protein is break across pbc. Using
        ``convert_trj()``
        with default parameters will fix it.

        :param nbox: list of 3 string for number of x, y, z dimensions copy
        :type nbox: list of string

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param cmd_args: Optional arguments for ``gmx genconf``

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.coor_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/copy_box/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../copy_box/top_SH3/1y0m_pdb2gmx.pdb -o \
.../copy_box/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        >>> prot.copy_box(nbox=[4,1,1])
        - Copy pbc box using genconf

        .. note::
            If ``name`` is not defined, the command will create a new pdb
            file name after the input one and adding "_copy_box.pdb".
            If ``name`` is defined the pdb filed will be saved in the same
            directory as input file, "_copy_box.pdb" will be added to
            ``name``.
        """

        logger.info("- Copy pbc box using genconf")
        # If name is not define use the object coor name and add _box.pdb
        if name is None:
            copy_coor = self.coor_file[:-4] + "_copy_box.pdb"
        # If not use the coor_path Path and add name and _box.pdb
        else:
            copy_coor = os.path.join(
                os_command.get_directory(self.coor_file),
                name + "_copy_box.pdb")
        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(copy_coor):
            logger.info("create_box not launched {} already exist".format(
                copy_coor))
            self.coor_file = copy_coor
            return

        if renumber:
            renum = 'yes'
        else:
            renum = 'no'

        # Nbox need to be a string list:
        nbox_str = [str(i) for i in nbox]

        cmd_copy = os_command.Command([GMX_BIN, "genconf",
                                       "-f", self.coor_file,
                                       "-o", copy_coor,
                                       "-renumber", renum,
                                       "-nbox"] + nbox_str,
                                      **cmd_args)

        cmd_copy.run()

        self.coor_file = copy_coor

    def solvate_box(self, out_folder, name=None, radius=0.21, cs=WATER_GRO,
                    check_file_out=True):
        """Solvate the pbc box with water or another mol defined with ``cs``
        using the ``gmx solvate`` command.


        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param radius: default van der Waals distance (nm)
        :type name: float, optional, default=0.21

        :param cs: solvant coordinate file
        :type cs: str, optional, default=``WATER_GRO``

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/solv_box/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../solv_box/top_SH3/1y0m_pdb2gmx.pdb -o \
.../solv_box/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        >>> prot.solvate_box(out_folder=TEST_OUT+'/solv_box/top_SH3_water/')
        - Solvate the pbc box
        Copy topologie file and dependancies

        .. note::
            If ``name`` is not defined, the command will create a new .pdb
            and .top file name after the object name and adding "_water".
        """

        logger.info("- Solvate the pbc box")

        # Create the out dir:
        start_dir = os.path.abspath(".")
        # Go in out_folder:
        os_command.create_and_go_dir(out_folder)

        if name is None:
            name = self.name + "_water"

        pdb_out = name + ".pdb"
        top_out = name + ".top"

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            logger.info("solvate_box not launched {} already exist".format(
                pdb_out))
            self.coor_file = pdb_out
            self.top_file = top_out
            os.chdir(start_dir)
            return

        # Copy the top file to the new directorie:
        topologie = TopSys(self.top_file)
        # topologie.display()
        topologie.copy_top_and_dependancies(top_out)

        cmd_solvate = os_command.Command([GMX_BIN, "solvate",
                                          "-cp", self.coor_file,
                                          "-cs", cs,
                                          "-p", top_out,
                                          "-o", pdb_out,
                                          "-radius", str(radius)])

        # cmd_solvate.display()
        cmd_solvate.run()

        self.coor_file = pdb_out
        self.top_file = top_out

        os.chdir(start_dir)

    def add_ions(self, out_folder, name=None, ion_C=0.15, pname="NA",
                 nname="CL", solv_name="SOL", maxwarn=1,
                 check_file_out=True):
        """Add ion in a system to neutralise the sys_charge and to reach the
        ionic concentration ``ion_C``.

        Ion number are computed using the water number and the charge of the
        system:

        1. With :math:`cation_{num} = {int(C_{ion} * water_{num}) \\over 55.5}`
        2. if :math:`cation_{num} + sys_{charge} >= 0` then\
            :math:`anion_{num} = cation_{num} + sys_{charge}` else \
            :math:`cation_{num} = -sys_{charge}`

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ion_C: ionic concentraton (Molar)
        :type ion_C: float, optional, default=0.15

        :param pname: cation name
        :type pname: str, optional, default="NA"

        :param pname: anion name
        :type pname: str, optional, default="CL"

        :param pname: solvant name
        :type pname: str, optional, default="SOL"

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/add_ions/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the \
topologie in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.create_box() #doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../add_ions/top_SH3/1y0m_pdb2gmx.pdb -o \
.../add_ions/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
        >>> prot.solvate_box(out_folder=TEST_OUT+'/add_ions/top_SH3_water/')
        - Solvate the pbc box
        Copy topologie file and dependancies
        >>> prot.add_ions(out_folder=TEST_OUT+'/add_ions/top_SH3_water_ions/')\
        #doctest: +ELLIPSIS
        Copy topologie file and dependancies
        - Create the tpr file genion_1y0m_ion.tpr
        gmx grompp -f .../template/mini.mdp -c \
../top_SH3_water/1y0m_water.pdb -r ../top_SH3_water/1y0m_water.pdb \
-p 1y0m_ion.top -po out_mini.mdp -o genion_1y0m_ion.tpr -maxwarn 1
        - Add ions to the system with an ionic concentration of 0.15 M , \
sytem charge = 0.0 water \
num= 56...
        Add ions : NA : 15   CL : 15
        gmx genion -s genion_1y0m_ion.tpr -p 1y0m_ion.top -o 1y0m_ion.gro \
-np 15 -pname NA -nn 15 \
-nname CL
        >>> prot.em(out_folder=TEST_OUT+'/add_ions/em_SH3_water_ions/', \
nsteps=10, constraints="none")
        - Create the tpr file 1y0m.tpr
        gmx grompp -f 1y0m.mdp -c ../top_SH3_water_ions/1y0m_ion.gro -r \
../top_SH3_water_ions/1y0m_ion.gro -p ../top_SH3_water_ions/1y0m_ion.top \
-po out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
        - Launch the simulation 1y0m.tpr
        gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright


        .. note::
            If ``name`` is not defined, the command will create a new .pdb
            and .top file name after the object name and adding "_ion".

        .. note ::
            There might be some charge issues with amber forcefield for
            example. There is a discrepency between the atom charge and
            the total charge column with amber. In the itp charge is chown
            with a 2 decimal precision as in the rtp file it can be up to
            5 decimals.
            Should consider using the total charge to deduce the atom charge
            and avoid errors. Up to now it has been fixed using `round`
            function instead of `int` for system charge
        """

        if name is None:
            name = self.name + "_ion"

        # Create the out dir:
        start_dir = os.path.abspath(".")
        # Go in out_folder:
        os_command.create_and_go_dir(out_folder)

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(
                name + ".gro"):
            logger.info("add ions not launched {} already exist".format(
                name + ".gro"))
            self.coor_file = name + ".gro"
            self.top_file = name + ".top"
            os.chdir(start_dir)
            return

        # Copy the top file to the new directorie:
        topologie = TopSys(self.top_file)
        topologie.copy_top_and_dependancies(name + ".top")

        self.top_file = name + ".top"

        # Create tpr:
        self.mdp = os.path.join(GROMACS_MOD_DIRNAME, "template/mini.mdp")
        self.add_tpr(name="genion_" + name, maxwarn=maxwarn)

        # Get charge:
        top = TopSys(self.top_file)
        sys_charge = top.charge()
        water_num = top.mol_num("SOL")
        logger.info("- Add ions to the system with an ionic concentration"
                    " of {} M , sytem charge = {} water num= {}".format(
                        ion_C, sys_charge, water_num))

        cation_num = round(ion_C / 55.5 * water_num)
        # Check if anion_num (cation_num  + sys_charge) is negative,
        # raise the canion_num to the sys_charge absolute value
        if (cation_num + sys_charge) < 0:
            cation_num = round(-1 * sys_charge)
        anion_num = round(cation_num + sys_charge)

        logger.info("Add ions : {} : {}   {} : {}".format(
            pname, cation_num, nname, anion_num))

        cmd_ions = os_command.Command([GMX_BIN, "genion",
                                       "-s", self.tpr,
                                       "-p", self.top_file,
                                       "-o", name + ".gro",
                                       "-np", str(cation_num),
                                       "-pname", pname,
                                       "-nn", str(anion_num),
                                       "-nname", nname])

        cmd_ions.display()
        cmd_ions.run(com_input=solv_name)

        self.coor_file = name + ".gro"

        os.chdir(start_dir)

    def solvate_add_ions(self, out_folder, name=None, ion_C=0.15,
                         create_box_flag=True, box_dist=1.1,
                         radius=0.25,
                         maxwarn=1):
        """Solvate a system with three succesive steps:

            1. Create box using ``create_box()``
            2. Add water using ``solvate_box()``
            3. Add ions using ``add_ions()``

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ion_C: ionic concentraton (Molar)
        :type ion_C: float, optional, default=0.15

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.coor_file solvate_box
            * self.top_file

        :Example:

        >>> TEST_OUT = getfixture('tmpdir')
        >>> prot = GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
        >>> prot.add_top(out_folder=TEST_OUT+'/solvate_add_ions/top_SH3/')\
        #doctest: +ELLIPSIS
        - Create topologie
        gmx pdb2gmx -f .../test_files/1y0m.pdb -o 1y0m_pdb2gmx.pdb -p \
1y0m_pdb2gmx.top -i 1y0m_posre.itp -water tip3p -ff charmm36-jul2017
        Molecule topologie present in 1y0m_pdb2gmx.top , extract the topologie\
 in a separate file: 1y0m_pdb2gmx.itp
        Protein_chain_A
        - ITP file: 1y0m_pdb2gmx.itp
        - molecules defined in the itp file:
        * Protein_chain_A
        Rewrite topologie: 1y0m_pdb2gmx.top
        >>> prot.solvate_add_ions(out_folder=TEST_OUT+'/solvate_add_ions/\
top_SH3_water_ions/')#doctest: +ELLIPSIS
        - Create pbc box
        gmx editconf -f .../solvate_add_ions/top_SH3/1y0m_pdb2gmx.pdb \
-o .../solvate_add_ions/top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron \
-d 1.1
        - Solvate the pbc box
        Copy topologie file and dependancies
        Copy topologie file and dependancies
        - Create the tpr file genion_1y0m_water_ion.tpr
        gmx grompp -f .../template/mini.mdp -c 1y0m_water.pdb -r \
1y0m_water.pdb -p 1y0m_water_ion.top -po out_mini.mdp -o \
genion_1y0m_water_ion.tpr -maxwarn 1
        - Add ions to the system with an ionic concentration of 0.15 M , \
sytem charge = 0.0 water num= 62...
        Add ions : NA : 17   CL : 17
        gmx genion -s genion_1y0m_water_ion.tpr -p 1y0m_water_ion.top -o \
1y0m_water_ion.gro -np 17 -pname NA -nn 17 -nname CL
        >>> prot.em(out_folder=TEST_OUT+'/solvate_add_ions/em_SH3_water_ions/'\
, nsteps=10, constraints = "none")
        - Create the tpr file 1y0m.tpr
        gmx grompp -f 1y0m.mdp -c ../top_SH3_water_ions/1y0m_water_ion.gro \
-r ../top_SH3_water_ions/1y0m_water_ion.gro -p \
../top_SH3_water_ions/1y0m_water_ion.top -po out_1y0m.mdp -o 1y0m.tpr \
-maxwarn 1
        - Launch the simulation 1y0m.tpr
        gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 \
-nocopyright



        .. note::
            If ``name`` is not defined, it will use the object name.


        """

        if name is None:
            name = self.name

        # Create box:
        if create_box_flag:
            self.create_box(dist=box_dist)

        # Solvate box:
        self.solvate_box(out_folder=out_folder,
                         name=name + "_water", radius=radius)

        # Add ions:
        self.add_ions(out_folder=out_folder,
                      name=name + "_water_ion",
                      ion_C=ion_C, maxwarn=maxwarn)

    def create_peptide(self, sequence, out_folder, N_ter="None", C_ter="COOH",
                       em_nsteps=1000, equi_nsteps=10000, posre_post="_pep",
                       vsite='none'):
        """Create a linear peptide structure and topologie:

            1. Create a peptide with pymol with one more residue G at the
                beginning of the peptide. This residue will then be change to
                an ACE. NH2 terminaison raise some issue with virtual sites
                and cannot be used.
            2. Create the topologie using ``add_top()``
            3. Minimise the structure using ``em()``
            4. Do a vacuum equilibration of the peptide using ``run_md_sim()``

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: generic name of the system
        :type name: str, optional, default=None

        :param ion_C: ionic concentraton (Molar)
        :type ion_C: float, optional, default=0.15

        :param vsite: option for topologie's bonds constraints ("none",
            "hydrogens", "all")
        :type vsite: str, optional, default="none"

        **Object requirement(s):**

            * None

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file


        :Example:

        .. code-block:: python

            > pep = GmxSys(name='SAM_pep')
            > pep.create_peptide(sequence='SAM', out_folder=os.path.join(\
str(TEST_OUT), 'peptide'), em_nsteps=10, equi_nsteps=10, vsite='hydrogens')
            -Make peptide: SAM
            residue name:X
            residue name:S
            residue name:A
            residue name:M
            Succeed to save file .../peptide/SAM.pdb
            - Create topologie
            gmx pdb2gmx -f ../SAM.pdb -o SAM_pdb2gmx.pdb -p SAM_pdb2gmx.top \
-i SAM_posre.itp -water tip3p -ff charmm36-jul2017 -ignh -ter -vsite hydrogens
            Molecule topologie present in SAM_pdb2gmx.top , extract the \
topologie in a separate file: SAM_pdb2gmx.itp
            Protein_chain_P
            - ITP file: SAM_pdb2gmx.itp
            - molecules defined in the itp file:
            * Protein_chain_P
            Rewrite topologie: SAM_pdb2gmx.top
            - Create pbc box
            gmx editconf -f .../peptide/00_top/SAM_pdb2gmx.pdb -o \
.../peptide/00_top/SAM_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
            - Create the tpr file SAM_pep.tpr
            gmx grompp -f SAM_pep.mdp -c ../00_top/SAM_pdb2gmx_box.pdb -r \
../00_top/SAM_pdb2gmx_box.pdb -p ../00_top/SAM_pdb2gmx.top -po \
out_SAM_pep.mdp -o SAM_pep.tpr -maxwarn 1
            - Launch the simulation SAM_pep.tpr
            gmx mdrun -s SAM_pep.tpr -deffnm SAM_pep -nt 0 -ntmpi 0 -nsteps \
-2 -nocopyright
            - Create the tpr file equi_vacuum_SAM.tpr
            gmx grompp -f equi_vacuum_SAM.mdp -c ../01_mini/SAM_pep.gro -r \
../01_mini/SAM_pep.gro -p ../00_top/SAM_pdb2gmx.top -po \
out_equi_vacuum_SAM.mdp -o equi_vacuum_SAM.tpr -maxwarn 1
            - Launch the simulation equi_vacuum_SAM.tpr
            gmx mdrun -s equi_vacuum_SAM.tpr -deffnm equi_vacuum_SAM -nt 0 \
-ntmpi 0 -nsteps -2 -nocopyright

        .. Warning::
            The peptide function won't work with gromacs version above 2018.
            There is issues with COOH C-temini, see:
            https://redmine.gromacs.org/issues/3301
            Use another C-ter or use a previous version of gromacs.


        """

        N_ter_dic = {"NH3+": "0", "NH2": "1", "5TER": "2", "None": "3"}
        C_ter_dic = {"COO-": "0", "COOH": "1", "CT2": "2", "3TER": "3",
                     "None": "4"}

        # Create a peptide with pymol with one more residue G at the beginning
        # of the peptide. This residue will then be change to an ACE
        # NH2 terminaison raise some issue with virtual sites and cannot be
        # used.
        pep_coor = pdb_manip.Coor()
        out_pdb = os.path.join(out_folder, sequence + ".pdb")
        pep_coor.make_peptide(sequence, out_pdb)
        self.coor_file = out_pdb

        self.add_top(out_folder=os.path.join(out_folder, "00_top"),
                     name=sequence, water="tip3p",
                     ff="charmm36-jul2017",
                     pdb2gmx_option_dict={'vsite': vsite,
                                          'ignh': None, 'ter': None},
                     check_file_out=False,
                     input_pdb2gmx=N_ter_dic[N_ter] + "\n" + C_ter_dic[C_ter],
                     posre_post=posre_post)

        # Create the box:
        # self.create_box()
        # Minimize the peptide:
        self.em(out_folder=os.path.join(out_folder, "01_mini"),
                nsteps=em_nsteps, constraints="none",
                create_box_flag=True)

        # Do sa short equi:
        if equi_nsteps > 0:

            if vsite != 'none':
                mdp_template = os.path.join(GROMACS_MOD_DIRNAME,
                                            "template/equi_vsites.mdp")
            else:
                mdp_template = os.path.join(GROMACS_MOD_DIRNAME,
                                            "template/equi.mdp")

            self.run_md_sim(out_folder=os.path.join(out_folder,
                                                    "02_equi_vacuum"),
                            name="equi_vacuum_" + sequence,
                            pdb_restr=None,
                            mdp_template=mdp_template,
                            maxwarn=1,
                            mdp_options={'nsteps': int(equi_nsteps),
                                         'dt': 0.001,
                                         'tc_grps': 'System',
                                         'tau_t': 0.1,
                                         'ref_t': 300,
                                         'pcoupl': 'no'})

    def insert_mol_sys(self, mol_gromacs, mol_num, new_name,
                       out_folder, check_file_out=True):
        """Insert a new molecule in a system:

        Insert structure and topologie of ``mol_num`` copy of ``mol_gromacs``
        molecule, in the system with 6 successive steps:

        1. Copy the molecule ``mol_num`` time.
        2. Change the chain ID of mol_gromacs to "Y", this step is necessary
            for vmd to recognize the inserted mol.
        3. Concat the two structure.
        4. Insert the molecule in the solvant with a vmd script.
        5. Update the topologie with the molecule and new water number.
        6. If the charge is not null add ions to neutralize the system.


        :param mol_gromacs: molecule object to be inserted
        :type mol_gromacs: GmxSys object

        :param mol_num: molecule number to be inserted
        :type mol_num: int

        :param new_name: generic name of the system
        :type new_name: str

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file
            * self.top_file

        **Object field(s) changed:**

            * self.coor_file
            * self.top_file

        .. note::
            VMD don't need anymore to be installed to run the peptide creation
        """

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(
                os.path.join(out_folder, new_name + ".top")):
            logger.info("insert_mol_sys not launched {} already exist".format(
                out_folder + "/" + new_name + ".top"))
            if os_command.check_file_and_create_path(
                    os.path.join(out_folder, new_name + "_neutral.pdb")):
                self.coor_file = os.path.join(
                    out_folder, new_name + "_neutral.pdb")
                self.top_file = os.path.join(out_folder, new_name + ".top")
                return
            elif os_command.check_file_and_create_path(
                    os.path.join(out_folder, new_name + ".pdb")):
                self.coor_file = os.path.join(out_folder, new_name + ".pdb")
                self.top_file = os.path.join(out_folder, new_name + ".top")
                return
            logger.error('Error top file exist but not coor file')
            raise IOError('coor file not found')

        # Create and got to the out dir:
        start_dir = os.path.abspath(".")
        os_command.create_and_go_dir(out_folder)

        # Get molecule resiude number:
        mol_gromacs.tpr = mol_gromacs.coor_file
        mol_gromacs.convert_trj(traj=False, pbc='none')
        mol_coor = pdb_manip.Coor(mol_gromacs.coor_file)
        res_num = len(mol_coor.get_attribute_selection(
            attribute='uniq_resid'))

        # Copy the mol using genconf:
        # Add random rotation ?
        if mol_num != 1:
            if res_num == 1:
                renum = True
            else:
                renum = False
            mol_gromacs.copy_box(nbox=[mol_num, 1, 1],
                                 check_file_out=check_file_out,
                                 rot="yes", renumber=renum)

        # Before doing the concat, Change the chain of mol_pdb to "Y",
        # this step is necessary for vmd to reognize the inserted mol
        mol_coor = pdb_manip.Coor(mol_gromacs.coor_file)
        mol_coor.change_pdb_field({"chain": "Y"})
        mol_coor.write_pdb(mol_gromacs.coor_file, check_file_out=False)
        # mol_length = int(mol_coor.get_aa_num() / mol_num)
        res_num = len(mol_coor.get_attribute_selection(
            attribute='uniq_resid'))
        logger.info("Res num: {}".format(res_num))

        # Concat the two pdb sys_pdb and mol_pdb
        concat_sys = new_name + "_pre_mix.pdb"

        # Get a compact pdb for the sys pdb, need to add a tpr if not already
        if self.tpr is None:
            self.sim_name = "tmp"
            mini_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                             "template/mini.mdp")
            self.add_mdp(mdp_template=mini_template_mdp, mdp_options={})
            self.add_tpr(name="tmp")
        self.convert_trj(traj=False)
        GmxSys.concat_coor(self.coor_file, mol_gromacs.coor_file,
                           pdb_out=concat_sys)

        # Do the molecule insertion with the pdb_manip module:

        sys_pdb = pdb_manip.Coor(concat_sys)

        sys_pdb.insert_mol(pdb_out=new_name + ".pdb", out_folder=".",
                           mol_chain="Y", mol_num=mol_num,
                           check_file_out=check_file_out)

        self.coor_file = new_name + ".pdb"

        # Insert the peptide top in the prot_sys top
        # Copy itp and posre files of mol_top to the new location
        top_mol = TopSys(mol_gromacs.top_file)
        old_name = top_mol.mol_comp[0]['name']
        # print("Old topologie name is:", old_name)
        top_mol.change_mol_name(old_name, "Ligand")
        top_mol.copy_dependancies("./")
        # top_mol.display()
        # Get the new location of the peptide itp file:
        mol_itp = top_mol.get_include_no_posre_file_list()
        # print("Include:", mol_itp)

        # Get the system topologie:
        sys_topologie = TopSys(self.top_file)
        # Add atomtypes:
        if mol_itp[0].endswith('atomtypes.itp'):
            sys_topologie.add_atomtypes(mol_itp[0])
            lig_itp = mol_itp[-1]
        else:
            lig_itp = mol_itp[0]
        # Add the peptide in the sys topologie and update the water num:
        sys_topologie.add_mol(mol_name="Ligand", mol_itp_file=lig_itp,
                              mol_num=mol_num)

        # Get the new water num after peptide insertion:
        sys_dict = pdb_manip.Coor(self.coor_file)
        water_res = sys_dict.get_attribute_selection(
            selec_dict={"res_name": ["SOL"]}, attribute='uniq_resid')
        logger.info("Water num: {}".format(len(water_res)))
        sys_topologie.change_mol_num(mol_name="SOL", mol_num=len(water_res))
        # save the top:
        sys_topologie.write_file(new_name + ".top")

        self.top_file = new_name + ".top"

        charge = sys_topologie.charge()
        logger.info("CHARGE: {}".format(charge))
        if charge != 0:
            if not os.path.isfile(new_name + "_neutral.pdb"):
                logger.info("Should neutralize the system")
                self.add_ions(out_folder=".", name=new_name + "_neutral",
                              ion_C=0)

        os.chdir(start_dir)

    def extract_mol_sys(self, out_folder, res_name):
        """Extract a molecule topologie and coordinates from a system:

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param res_name: Molecule residue name to extract
        :type res_name: str

        :return: The GmxSys of the molecule alone
        :rtype: GmxSys

        .. note::
            This function does not use or affect the GmxSys object.

        """

        # Extract coordinates:
        if self.coor_file[-3:] != 'pdb':
            self.convert_trj(traj=False)
        full_coor = pdb_manip.Coor(self.coor_file)
        # Select coor:
        mol_coor = full_coor.select_part_dict(
            selec_dict={'res_name': [res_name]})
        mol_coor.write_pdb(os.path.join(out_folder, res_name+'_only.pdb'))

        # Extract topologie:

        # Get the system topologie:
        sys_topologie = TopSys(self.top_file)
        mol_top = copy.deepcopy(sys_topologie)

        itp_list = []
        # Keep only res_name itp
        for itp in mol_top.itp_list:

            # print(itp.name)
            to_keep = False
            mol_num = 0

            for top_mol in itp.top_mol_list:
                mol_num += 1
                # print('\t'+top_mol.name)
                if top_mol.name in [res_name, 'SOL']:
                    to_keep = True
                    break
            if to_keep or mol_num == 0:
                itp_list.append(itp)

        mol_top.itp_list = itp_list

        new_mol_comp = []
        # Keep only res_name in mol composition
        for mol in mol_top.mol_comp:
            if mol['name'] == res_name:
                new_mol_comp.append(mol)

        mol_top.mol_comp = new_mol_comp

        mol_top.display()
        mol_top.write_file(os.path.join(out_folder, res_name+'_only.top'))

        mol_sys = GmxSys(name=res_name,
                         coor_file=os.path.join(out_folder,
                                                res_name+'_only.pdb'),
                         top_file=os.path.join(out_folder,
                                               res_name+'_only.top'))
        return mol_sys

    @staticmethod
    def concat_coor(*coor_in_files, pdb_out, check_file_out=True):
        """Concat a list of coordinates file in one coordinate file:

        :param coor_in_files: list of pdb/gro files
        :type coor_in_files: list of str

        :param pdb_out: file to save the concat structure
        :type pdb_out: str

        :return: name of the new pdb file
        :rtype: str

        .. note::
            This function does not use or affect the GmxSys object.

        """

        # Check if output files exist:
        if check_file_out and os.path.isfile(pdb_out):
            logger.info("PDB files not created, {} already exist".format(
                pdb_out))
            return

        pdb_in_files = []

        for coor_in in coor_in_files:
            # print("File:", coor_in)
            if (coor_in[-3:]) == "pdb":
                pdb_in_files.append(coor_in)
            elif (coor_in[-3:]) == "gro":
                tmp_gromacs = GmxSys(coor_file=coor_in, tpr=coor_in)
                tmp_gromacs.convert_trj(traj=False, pbc='none')
                pdb_in_files.append(tmp_gromacs.coor_file)
            else:
                raise RuntimeError('Cannot concat the file, should be'
                                   ' gro or pdb format')
        logger.info("Concat files: {}".format(pdb_in_files))
        return pdb_manip.Coor.concat_pdb(pdb_out=pdb_out, *pdb_in_files)

    def concat_traj(self, *xtc_files_list, concat_traj_out,
                    check_file_out=True):
        """Concat a list of trajectory file in one trajectory file:

        :param xtc_in_files: list of xtc files
        :type xtc_in_files: list of str

        :param concat_traj_out: file to save the concat trajectory
        :type concat_traj_out: str


        **Object field(s) changed:**

            * self.xtc
        """

        # Check if output files exist:
        if check_file_out and os.path.isfile(concat_traj_out):
            logger.info("XTC files not created, {} already exist".format(
                concat_traj_out))
            self.xtc = concat_traj_out
            return

        cmd_list = [GMX_BIN, "trjcat",
                    "-o", concat_traj_out,
                    "-f"] + list(xtc_files_list)

        cmd_trjcat = os_command.Command(cmd_list)
        cmd_trjcat.display()
        cmd_trjcat.run(display=False)

        self.xtc = concat_traj_out

        return

    def concat_edr(self, *edr_files_list, concat_edr_out, check_file_out=True):
        """Concat a list of energy file in one energy file:

        :param xtc_in_files: list of edr files
        :type xtc_in_files: list of str

        :param concat_edr_out: file to save the concat energy
        :type concat_edr_out: str


        **Object field(s) changed:**

            * self.edr
        """
        # Check if output files exist:
        if check_file_out and os.path.isfile(concat_edr_out):
            logger.info("Edr files not created, {} already exist".format(
                concat_edr_out))
            self.edr = concat_edr_out
            return

        cmd_list = [GMX_BIN, "eneconv",
                    "-o", concat_edr_out,
                    "-f"] + list(edr_files_list)

        cmd_eneconv = os_command.Command(cmd_list)
        cmd_eneconv.display()
        cmd_eneconv.run(display=False)

        self.edr = concat_edr_out

        return

    ##########################################################
    # ###########  SIMULATION RELATED FUNCTIONS  #############
    ##########################################################

    def get_mdp_dict(self, mdp_file=None):
        """ Extract mdp file's parameters and return it as a dict.
        By default read self.coor, but if mdp_file option is defined,
        it will read it.

        :param mdp_file: mdp file
        :type mdp_file: str, default=None

        **Object requirement(s):**

            * self.mdp

        """

        if mdp_file is None:
            mdp_file = self.mdp

        mdp_dict = {}

        with open(mdp_file) as filein:

            for line in filein:
                if line[0] != ';':
                    line_no_space = line[:-1].replace(" ", "")
                    line_split = line_no_space.split('=')
                    # print(line)
                    # print(line_split)
                    if len(line_split) == 2:
                        key = line_split[0].lower().replace("-", "_")
                        mdp_dict[key] = line_split[1].lower()

        return(mdp_dict)

    def add_mdp(self, mdp_template, mdp_options, folder_out="",
                check_file_out=True):
        """Create the MD simulation input mdp file from a template.

        Read a template mdp file and replace define fields in mdp_options
        with the new value. In case the field name has a '-' , repalce it
        by : '_'.

        :param mdp_template: mdp file template
        :type mdp_template: str

        :param mdp_options: New parameters to use
        :type mdp_options: dict

        :param folder_out: Path for output file
        :type folder_out: str, default=""

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.sim_name

        **Object field(s) changed:**

            * self.mdp

        """

        mdp_out = self.sim_name + ".mdp"

        if folder_out != "":
            mdp_out = os.path.join(folder_out, mdp_out)

        # Check if output files exist:
        if check_file_out and os.path.isfile(mdp_out):
            logger.info("Mdp files not created, {} already exist".format(
                mdp_out))
            self.mdp = mdp_out
            return

        # Replace - by _:
        local_mdp_opt = {}
        for key, value in mdp_options.items():
            local_mdp_opt[key.replace("-", "_").lower()] = str(value)

        # local_mdp_opt = mdp_options.copy()

        mdp_template_dict = self.get_mdp_dict(mdp_template)

        # Check options not founded is the mdp_template
        for key, value in local_mdp_opt.items():
            if key not in mdp_template_dict:
                logger.warning("WARNING !!! ADDING unusual parameter : \"{}"
                               "\"in the mdp file {}".format(key, mdp_out))

        mdp_template_dict.update(local_mdp_opt)

        # Write mdp file:
        filout = open(mdp_out, 'w')
        for key, value in mdp_template_dict.items():
            line = "  {:25}   =    {}\n".format(
                key.replace("-", "_").lower(), str(value))
            filout.write(line)

        filout.close()

        self.mdp = mdp_out

    def create_mdp(self, mdp_options, folder_out="", check_file_out=True):
        """Create the MD simulation input mdp file.

        :param mdp_options: New parameters to use
        :type mdp_options: dict

        :param folder_out: Path for output file
        :type folder_out: str, default=""

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.sim_name

        **Object field(s) changed:**

            * self.mdp

        """

        mdp_out = self.sim_name + ".mdp"

        if folder_out != "":
            mdp_out = os.path.join(folder_out, mdp_out)

        # Check if output files exist:
        if check_file_out and os.path.isfile(mdp_out):
            logger.info("Mdp files not created, {} already exist".format(
                mdp_out))
            self.mdp = mdp_out
            return

        filout = open(mdp_out, 'w')

        for key, value in mdp_options.items():
            line = "  {:25}   =    {}\n".format(
                key.replace("-", "_").lower(), str(value))
            filout.write(line)

        filout.close()

        self.mdp = mdp_out

    def add_ndx(self, ndx_cmd_input, ndx_name=None, folder_out="",
                check_file_out=True):
        """Create a ndx file using ``gmx make_ndx``

        :param ndx_cmd_input: Input arguments for ``gmx make_ndx``
        :type ndx_cmd_input: str

        :param ndx_name: output name for the index file
        :type ndx_name: str, default=None

        :param folder_out: Path for output file
        :type folder_out: str, default=""

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.coor_file

        **Object field(s) changed:**

            * self.ndx

        .. note::
            If name is not defined, will use the name of the object.
        """

        if ndx_name is not None:
            ndx_out = ndx_name + ".ndx"
        else:
            if self.tpr is not None:
                ndx_out = self.tpr[:-4] + ".ndx"
            else:
                ndx_out = self.coor_file[:-4] + ".ndx"

        if folder_out != "":
            ndx_out = os.path.join(folder_out, ndx_out)

        logger.info("- Create the ndx file {}".format(ndx_out))

        # Check if output files exist:
        if check_file_out and os.path.isfile(ndx_out):
            logger.info("add_ndx not launched {} already exist".format(
                ndx_out))
            self.ndx = ndx_out
            return

        # Add NDX input ?
        #             "-n", self.ndx,

        cmd_list = [GMX_BIN, "make_ndx",
                    "-f", self.coor_file,
                    "-o", ndx_out]

        cmd_ndx = os_command.Command(cmd_list)

        cmd_ndx.display()

        cmd_ndx.run(com_input=ndx_cmd_input, display=False)

        self.ndx = ndx_out

        # Clean index and remove multiple entry for the same
        # groups:
        with open(self.ndx, 'r') as ndx_file:
            line_list = ndx_file.readlines()

        group_list = []
        with open(self.ndx, 'w') as ndx_file:
            for line in line_list:
                if line.startswith('['):
                    group = line.strip()[1:-1].strip()
                    if group in group_list:
                        save_group = False
                    else:
                        save_group = True
                        group_list.append(group)
                if save_group:
                    ndx_file.write(line)

    def get_index_dict(self):
        """ Read and `.ndx` file and return and dictionnary
        with keys being the group name, nad value the index.
        """

        # Add ndx file if absent:
        if self.ndx is None:
            self.add_ndx('q\n')

        index_dict = {}
        index = 0
        with open(self.ndx, 'r') as ndx_file:
            line_list = ndx_file.readlines()
            for line in line_list:
                if line.startswith('['):
                    group = line.strip()[1:-1].strip()
                    index_dict[group] = index
                    index += 1
        return(index_dict)

    def center_mol_box(self, sele_dict={'res_name': PROT_RES + DNARNA_RES,
                                        'name': HA_NAME},
                       traj=False, ref_coor=None, **cmd_args):
        """ Center a sytem on a selection of residue

        :param res_list: List of residues
        :type res_list: str

        """

        if ref_coor is None:
            ref_coor = self.coor_file

        coor_complex = pdb_manip.Coor(ref_coor)
        center_res = coor_complex.get_center_residue(
            selec_dict=sele_dict,
            field='uniq_resid')
        # Gromacs residue index starts at 1
        center_res += 1

        self.add_ndx('ri {} \n q \n'.format(center_res), check_file_out=False)

        index_dict = self.get_index_dict()
        sel_center = index_dict['r_{}'.format(center_res)]

        self.convert_trj(select='{:d} \n System'.format(sel_center), traj=traj,
                         ur='tric',
                         pbc='mol', center='yes', **cmd_args)
        # self.convert_trj(traj=traj)

    def add_tpr(self, name, r=None, po=None, folder_out="",
                check_file_out=True, **grompp_options):
        r"""Create a tpr file using ``gmx grompp``.

        :param name: name of the simulation
        :type name: str

        :param r: reference coordinate file for position restraints
        :type r: str, default=None

        :param po: output file for the mdp file
        :type po: str, default=None

        :param folder_out: Path for output file
        :type folder_out: str, default=""

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param \**grompp_options: Optional arguments for ``gmx grompp``

        **Object requirement(s):**

            * self.mdp
            * self.coor_file
            * self.top_file

        **Object optional requirement(s):**

            * self.ndx

        **Object field(s) changed:**

            * self.tpr

        """

        self.sim_name = name
        tpr_out = self.sim_name + ".tpr"

        if folder_out != "":
            tpr_out = os.path.join(folder_out, tpr_out)

        logger.info("- Create the tpr file {}".format(tpr_out))

        # Check if output files exist:
        if check_file_out and os.path.isfile(tpr_out):
            logger.info("create_tpr not launched {} already exist".format(
                tpr_out))
            self.tpr = tpr_out
            return

        # Define grompp few input variable as function of previous input:
        if r is None:
            r = self.coor_file

        if po is None:
            po = "./out_" + self.mdp.split("/")[-1]

        cmd_list = [GMX_BIN, "grompp",
                    "-f", self.mdp,
                    "-c", self.coor_file,
                    "-r", r,
                    "-p", self.top_file,
                    "-po", po,
                    "-o", tpr_out]

        for key, value in grompp_options.items():
            cmd_list = cmd_list + ["-" + key, str(value)]

        # Add ndx if defined in the object
        if self.ndx is not None:
            cmd_list = cmd_list + ["-n", self.ndx]

        cmd_tpr = os_command.Command(cmd_list)

        cmd_tpr.display()
        # cmd_tpr.define_env(my_env={**os.environ, 'GMXLIB': FORCEFIELD_PATH})
        cmd_tpr.define_env(my_env=os.environ.update({'GMXLIB': GMXLIB_var}))

        cmd_tpr.run(display=False)

        self.tpr = tpr_out

    def run_simulation(self, check_file_out=True, cpi=None, nsteps=-2,
                       rerun=False, monitor_tool=monitor.PROGRESS_BAR):
        """
        Launch the simulation using ``gmx mdrun``

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :param cpi: checkpoint file, if defined, it will restart a simulation
            and run ``nsteps``.
        :type cpi: str, default=None

        :param nsteps: Number of steps to run, (-2 : will use mdp parameter)
        :type nsteps: int, default=-2

        :param rerun: option to rerun a simulation (eg. recompute energy)
        :type rerun: bool, default=False

        :param monitor_tool: option to monitor a simulation, if not `None`,
            monitor should contains two values: ``function`` the function
            to be ran while simulation is running and ``input`` parameters
            for the function.
        :type monitor_tool: dict, default=None

        **Object requirement(s):**

            * self.tpr
            * self.sim_name
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Additional requirement(s) for rerun:**

            * self.xtc

        **Object field(s) changed:**

            * self.coor_file
            * self.xtc
            * self.edr
            * self.log

        .. note::
            The function must be launched in the path where the tpr is present.

        .. note::
            If cpi file is defined the simulation will restart with the
            ``-noappend`` option, if cpi is not defined, but the .cpt file
            exist, it will restart with "append".
        """

        # nsteps = -2 , will use the mdp file option

        logger.info("- Launch the simulation {}".format(self.tpr))

        # Check if output files exist:
        if check_file_out and os.path.isfile(self.sim_name + ".gro"):
            logger.info("Simulation not launched {} "
                        "already exist".format(self.sim_name + ".gro"))
            self.coor_file = self.sim_name + ".gro"
            if os_command.check_file_exist(self.sim_name + ".xtc"):
                self.xtc = self.sim_name + ".xtc"
            else:
                self.xtc = self.sim_name + ".trr"
            self.edr = self.sim_name + ".edr"
            self.log = self.sim_name + ".log"
            if os_command.check_file_exist(self.sim_name + ".xvg"):
                self.xvg = self.sim_name + ".xvg"
            return

        if rerun and check_file_out and os.path.isfile(
                self.sim_name + ".edr"):
            logger.info("Simulation not launched {} "
                        "already exist".format(self.sim_name + ".edr"))
            self.edr = self.sim_name + ".edr"
            return

        # Removing the "-v" option, seems to be harmless, I keep it for
        # subproccess reason
        # I prefer to get error with short test because of std output size
        # than with very long run after a long time and not seeing an error
        # is going one
        cmd_list = [GMX_BIN, "mdrun",
                    "-s", self.tpr,
                    "-deffnm", self.sim_name,
                    "-nt", str(int(self.nt)),
                    "-ntmpi", str(int(self.ntmpi)),
                    "-nsteps", str(int(nsteps)),
                    "-nocopyright"]

        # Don't add gpu_id if not specified
        if self.gpu_id is not None:
            cmd_list = cmd_list + ["-gpu_id", self.gpu_id]
        if rerun:
            cmd_list = cmd_list + ["-rerun", self.xtc]

        # If cpi file is included, do a restart from the cpi and don't append
        # xtc, log, edr, ...
        # If cpi option is not defined and the cpt exist, do a restart with
        # append
        if cpi is not None:
            cmd_list = cmd_list + ["-noappend"]
            cmd_list = cmd_list + ["-cpi", cpi]
        elif os.path.isfile(self.sim_name + ".cpt"):
            cmd_list = cmd_list + ["-append"]
            cmd_list = cmd_list + ["-cpi", self.sim_name + ".cpt"]

        cmd_run = os_command.Command(cmd_list)

        cmd_run.display()
        if monitor_tool is None:
            cmd_run.run()
        else:
            monitor_files = {'xtc': self.sim_name + ".xtc",
                             'edr': self.sim_name + ".edr",
                             'log': self.sim_name + ".log",
                             'nsteps': nsteps}
            sim_mdp_dict = self.get_mdp_dict()
            if nsteps == -2:
                monitor_files['nsteps'] = int(sim_mdp_dict['nsteps'])
            # Check if it is an extend simulation:
            # And need to update the log file to XXX.part0XXX
            elif os.path.isfile(self.sim_name + ".cpt"):
                monitor_files['nsteps'] += self.get_simulation_time() / float(
                    sim_mdp_dict['dt'])
                # Find last job number (sim_num)
                for sim_num in range(2, 9999):
                    if not os.path.isfile(
                            f"{self.sim_name}.part{sim_num:04d}.edr"):
                        break
                monitor_files['xtc'] = f"{self.sim_name}.part{sim_num:04d}.xtc"
                monitor_files['edr'] = f"{self.sim_name}.part{sim_num:04d}.edr"
                monitor_files['log'] = f"{self.sim_name}.part{sim_num:04d}.log"

            monitor_tool.update(monitor_files)
            cmd_run.run_background(monitor_tool)

        # If it's not a rerun, assign all output to the object variables
        # xtc, edr, log
        self.edr = self.sim_name + ".edr"
        if not rerun:
            self.coor_file = self.sim_name + ".gro"
            if os_command.check_file_exist(self.sim_name + ".xtc"):
                self.xtc = self.sim_name + ".xtc"
            else:
                self.xtc = self.sim_name + ".trr"
            self.log = self.sim_name + ".log"
            if os_command.check_file_exist(self.sim_name + ".xvg"):
                self.xvg = self.sim_name + ".xvg"

    def run_md_sim(self, out_folder, name, mdp_template, mdp_options,
                   pdb_restr=None, maxwarn=0,
                   monitor_tool=monitor.PROGRESS_BAR):
        """Run a simulation using 3 steps:

        1. Create a mdp file
        2. Create a tpr file using ``gmx grompp``
        3. Launch the simulation using ``gmx mdrun``

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str

        :param mdp_template: mdp file template
        :type mdp_template: str

        :param mdp_options: New parameters to use
        :type mdp_options: dict

        :param pdb_restr: reference coordinate file for position restraints
        :type pdb_restr: str, default=None

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id


        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc

        """

        start_dir = os.path.abspath(".")

        if pdb_restr is not None:
            pdb_restr = os_command.full_path_and_check(pdb_restr)

        if self.top_file is None:
            raise ValueError("Simulation could not be ran because the"
                             " topologie is not defined")

        # Create and go in out_folder:
        os_command.create_and_go_dir(out_folder)

        # Save previous state:
        self.save_state()

        # Create mdp :
        self.sim_name = name
        self.add_mdp(mdp_template=mdp_template, mdp_options=mdp_options)
        self.add_tpr(name=name, r=pdb_restr, maxwarn=maxwarn)
        self.run_simulation(monitor_tool=monitor_tool)

        # Get absolute path:
        os.chdir(start_dir)

    def em(self, out_folder, name=None, posres="",
           create_box_flag=False, monitor_tool=monitor.PROGRESS_BAR, maxwarn=1,
           **mdp_options):
        """Minimize a system.

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param nsteps: number of minimisation steps
        :type nsteps: int, default=1000

        :param posres: option for the ``define`` variable in the mdp file,
            need to be define to have postion restraints.
        :type posres: str, default=""

        :param create_box_flag: flag to create or not a box to the input coor
            file.
        :type create_box_flagt: bool, optional, default=False

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if create_box_flag:
            self.create_box()
        if name is None:
            name = self.name
        # mdp :
        mini_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                         "template/mini.mdp")
        mdp_options.update({'define': posres})

        self.run_md_sim(out_folder=out_folder, name=name,
                        mdp_template=mini_template_mdp,
                        monitor_tool=monitor_tool, mdp_options=mdp_options,
                        maxwarn=maxwarn)

    def em_2_steps(self, out_folder, name=None, no_constr_nsteps=1000,
                   constr_nsteps=1000,
                   posres="", create_box_flag=False,
                   monitor_tool=monitor.PROGRESS_BAR,
                   maxwarn=1,
                   **mdp_options):
        """Minimize a system in two steps:

        1. minimisation without bond constraints
        2. minimisation using bond constraint for bonds involving hydrogen

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param no_constr_nsteps: number of minimisation steps in the first step
        :type no_constr_nsteps: int, default=1000

        :param constr_nsteps: number of minimisation steps in the second step
        :type constr_nsteps: int, default=1000

        :param posres: option for the ``define`` variable in the mdp file,
            need to be define to have postion restraints.
        :type posres: str, default=""

        :param create_box_flag: flag to create or not a box to the input coor
            file.
        :type create_box_flag: bool, optional, default=False

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if name is None:
            name = self.name

        self.em(out_folder=out_folder, name="Init_em_" + name,
                nsteps=int(no_constr_nsteps),
                posres=posres, create_box_flag=create_box_flag,
                constraints="none",
                monitor_tool=monitor_tool,
                maxwarn=maxwarn,
                **mdp_options)

        self.em(out_folder=out_folder, name=name, nsteps=int(constr_nsteps),
                posres=posres, create_box_flag=False, constraints="all-bonds",
                maxwarn=maxwarn,
                monitor_tool=monitor_tool, **mdp_options)

    def equi_three_step(self, out_folder, name=None, pdb_restr=None,
                        nsteps_HA=100000,
                        nsteps_CA=200000, nsteps_CA_LOW=400000, dt=0.002,
                        dt_HA=0.001,
                        maxwarn=0,
                        monitor_tool=monitor.PROGRESS_BAR,
                        vsite='none', **mdp_options):
        """Equilibrate a system in 3 steps:

        1. equilibration of nsteps_HA with position restraints on Heavy Atoms
        with dt = dt_HA
        2. equilibration of nsteps_CA with position restraints on Carbon Alpha
        with dt = dt
        3. equilibration of nsteps_CA_LOW with position restraints on Carbon
        Alpha with Low restraints with dt = dt

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param pdb_restr: reference coordinate file for position restraints
        :type pdb_restr: str, default=None

        :param nsteps_HA: number of equilibration steps with HA constraints
        :type nsteps_HA: int, default=100000

        :param nsteps_CA: number of equilibration steps with CA constraints
        :type nsteps_CA: int, default=200000

        :param nsteps_CA_LOW: number of equilibration steps with CA_LOW
            constraints
        :type nsteps_CA_LOW: int, default=400000

        :param dt_HA: integration time step for HA equilibration
        :type dt_HA: float, default=0.001

        :param dt: integration time step for CA and CA_LOW equilibration
        :type dt: float, default=0.002

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param vsite: option for bonds constraints ("none")
        :type vsite: str, optional, default="none"

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc

        .. Note::
            In case of LINCS warning or segmentation fault, try to center the
            protein in the box using the center_mol_box() function.

        """

        if name is None:
            name = self.name

        if pdb_restr is None:
            pdb_restr = self.coor_file

        if vsite != 'none':
            equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                             "template/equi_vsites.mdp")
        else:
            equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                             "template/equi.mdp")
            if dt > 0.002 or dt_HA > 0.002:
                logger.warning('Be carfull, dt = {}/{} witout using vsite'
                               ', that may induce sever issues with your '
                               'simualtion !'.format(dt, dt_HA))

        # need to center prot/DNA in the box to avoid position restraint
        # across the box
        # self.center_mol_box()

        mdp_options.update({'nsteps': int(nsteps_HA),
                            'pcoupl': 'berendsen',
                            'define': '-DPOSRES',
                            'dt': dt_HA,
                            'comm_mode': 'none',
                            'refcoord_scaling': 'com'})

        self.run_md_sim(out_folder=os.path.join(out_folder, "00_equi_HA"),
                        name="equi_HA_" + name,
                        pdb_restr=pdb_restr, mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

        mdp_options.update({'nsteps': int(nsteps_CA),
                            'define': '-DPOSRES_CA', 'dt': dt})
        self.run_md_sim(out_folder=os.path.join(out_folder, "01_equi_CA"),
                        name="equi_CA_" + name,
                        pdb_restr=pdb_restr, mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

        mdp_options.update({'nsteps': int(nsteps_CA_LOW),
                            'pcoupl': 'parrinello-rahman',
                            'define': '-DPOSRES_CA_LOW'})
        self.run_md_sim(out_folder=os.path.join(out_folder, "02_equi_CA_LOW"),
                        name="equi_CA_LOW_" + name,
                        pdb_restr=pdb_restr, mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

    def em_equi_three_step_iter_error(self, out_folder, name=None,
                                      no_constr_nsteps=1000,
                                      constr_nsteps=1000,
                                      pdb_restr=None, nsteps_HA=100000,
                                      nsteps_CA=200000, nsteps_CA_LOW=400000,
                                      dt=0.002, dt_HA=0.001, maxwarn=0,
                                      iter_num=3,
                                      monitor_tool=monitor.PROGRESS_BAR,
                                      vsite='none',
                                      **mdp_options):
        """ Minimize a system in 2 steps:

        1. minimisation without bond constraints
        2. minimisation using bond constraint for bonds involving hydrogen

        Equilibrate a system in 3 steps:

        1. equilibration of nsteps_HA with position restraints on Heavy Atoms
        with dt = dt_HA
        2. equilibration of nsteps_CA with position restraints on Carbon Alpha
        with dt = dt
        3. equilibration of nsteps_CA_LOW with position restraints on Carbon
        Alpha with Low restraints with dt = dt

        In case this process will crash (eg. LINCS WARNING ...),
        the process will be rerun for `iter_num` time.

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param no_constr_nsteps: number of minimisation steps in the first step
        :type no_constr_nsteps: int, default=1000

        :param constr_nsteps: number of minimisation steps in the second step
        :type constr_nsteps: int, default=1000

        :param pdb_restr: reference coordinate file for position restraints
        :type pdb_restr: str, default=None

        :param nsteps_HA: number of equilibration steps with HA constraints
        :type nsteps_HA: int, default=100000

        :param nsteps_CA: number of equilibration steps with CA constraints
        :type nsteps_CA: int, default=200000

        :param nsteps_CA_LOW: number of equilibration steps with CA_LOW
            constraints
        :type nsteps_CA_LOW: int, default=400000

        :param dt_HA: integration time step for HA equilibration
        :type dt_HA: float, default=0.001

        :param dt: integration time step for CA and CA_LOW equilibration
        :type dt: float, default=0.002

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param vsite: option for bonds constraints ("none")
        :type vsite: str, optional, default="none"

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        start_sys = copy.deepcopy(self)
        start_dir = os.path.abspath(".")

        if vsite == 'none' and dt > 0.002 or dt_HA > 0.002:
            logger.warning('Be carfull, dt = {}/{} witout using vsite'
                           ', that may induce sever issues with your '
                           'simualtion !'.format(dt, dt_HA))

        for iter in range(iter_num):
            try:
                local_out_folder = out_folder + "/sys_em/"
                self.em_2_steps(local_out_folder, name=name,
                                no_constr_nsteps=no_constr_nsteps,
                                constr_nsteps=constr_nsteps,
                                posres="", create_box_flag=False,
                                monitor_tool=monitor_tool,
                                maxwarn=maxwarn,
                                **mdp_options)
                self.convert_trj(traj=False)

                local_out_folder = out_folder + "/sys_equi/"
                self.equi_three_step(local_out_folder, name=name,
                                     pdb_restr=pdb_restr, nsteps_HA=nsteps_HA,
                                     nsteps_CA=nsteps_CA,
                                     nsteps_CA_LOW=nsteps_CA_LOW,
                                     dt=dt, dt_HA=dt_HA,
                                     maxwarn=maxwarn,
                                     monitor_tool=monitor_tool,
                                     vsite=vsite, **mdp_options)
                break

            except RuntimeError as e:
                logger.error('Run {}/{} failed because of: {}'.format(
                    iter + 1, iter_num, e.args[0]))
                os.chdir(start_dir)
                # Remove directories
                for dir_to_del in [out_folder + "/sys_em/",
                                   out_folder + "/sys_equi/"]:
                    if os_command.check_directory_exist(dir_to_del):
                        os_command.delete_directory(dir_to_del)
                self = copy.deepcopy(start_sys)

        os.chdir(start_dir)

    def production(self, out_folder, name=None, nsteps=400000, dt=0.002,
                   maxwarn=0,
                   monitor_tool=monitor.PROGRESS_BAR, vsite='none',
                   **mdp_options):
        """Run a production run.

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param nsteps: number of minimisation steps
        :type nsteps: int, default=400000

        :param dt: number of minimisation steps
        :type dt: float, default=0.002

        :param maxwarn: Maximum number of warnings when using ``gmx grompp``
        :type maxwarn: int, default=0

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param vsite: option for topologie's bonds constraints ("none",
            "hydrogens", "all")
        :type vsite: str, optional, default="none"


        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if name is None:
            name = self.name

        if vsite != 'none':
            equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                             "template/equi_vsites.mdp")
        else:
            equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                             "template/equi.mdp")

        mdp_options.update({'nsteps': int(nsteps), 'dt': dt, 'define': ''})
        self.run_md_sim(out_folder=out_folder, mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, name="prod_" + name,
                        monitor_tool=monitor_tool, maxwarn=maxwarn)

    def extend_sim(self, tpr_file=None, nsteps=200000,
                   monitor_tool=monitor.PROGRESS_BAR):
        """Extend a simulation run.

        :param tpr_file: path of the tpr file
        :type tpr_path: str

        :param nsteps: number of steps
        :type nsteps: int, default=200000

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        **Object requirement(s):**

            * self.tpr
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.tpr
            * self.sim_name
            * self.coor_file
            * self.xtc

        .. warning::
            The simulation can only run nsteps more steps.
            Should be improved to check how many steps have been run, \
            and give a total number of step to compute.
        """

        if tpr_file is not None:
            self.tpr = tpr_file
        else:
            if self.tpr is None:
                logger.info("self.tpr or tpr_file is not defined \n"
                            "Could not restart the simulation")
                raise ValueError()

        # Get simulation time :
        sim_time = self.get_simulation_time()
        dt = float(self.get_mdp_dict()['dt'])
        nsteps_to_run = int(nsteps - sim_time / dt)
        if nsteps_to_run <= 0:
            logger.info("Simulation {} has already run"
                        " {} ps, extending simulation is useless.".format(
                            self.tpr[:-4], sim_time))
            return
        logger.info("- Extend simulation for {} steps".format(nsteps_to_run))

        self.sim_name = self.tpr.split("/")[-1][:-4]
        out_folder = os_command.get_directory(self.tpr)

        start_dir = os.path.abspath(".")

        os_command.create_and_go_dir(out_folder)

        self.run_simulation(cpi=self.sim_name + ".cpt",
                            nsteps=int(nsteps_to_run),
                            check_file_out=False, monitor_tool=monitor_tool)
        self.get_last_output()

        os.chdir(start_dir)

    def em_CG(self, out_folder, name=None, nsteps=500000,
              maxwarn=0,
              monitor_tool=monitor.PROGRESS_BAR,
              **mdp_options):
        """Equilibrate a system a CG system:

        1. equilibration of nsteps_HA with position restraints on Heavy Atoms
        with dt = dt_HA
        2. equilibration of nsteps_CA with position restraints on Carbon Alpha
        with dt = dt
        3. equilibration of nsteps_CA_LOW with position restraints on Carbon
        Alpha with Low restraints with dt = dt

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param pdb_restr: reference coordinate file for position restraints
        :type pdb_restr: str, default=None

        :param nsteps: number of equilibration steps with BB constraints
        :type nsteps: int, default=100000

        :param dt: integration time step for BB equilibration
        :type dt: float, default=0.002

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if name is None:
            name = self.name

        equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                         "template/martini_v2.x_new-rf.mdp")

        mdp_options.update({'integrator': 'steep', 'tcoupl': 'no',
                            'nsteps': int(nsteps)})

        self.run_md_sim(out_folder=out_folder, name="em_CG_" + name,
                        mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

    def equi_CG(self, out_folder, name=None, pdb_restr=None, nsteps=500000,
                dt=0.02, maxwarn=0,
                monitor_tool=monitor.PROGRESS_BAR,
                **mdp_options):
        """Equilibrate a system a CG system:

        1. equilibration of nsteps_HA with position restraints on Heavy Atoms
        with dt = dt_HA
        2. equilibration of nsteps_CA with position restraints on Carbon Alpha
        with dt = dt
        3. equilibration of nsteps_CA_LOW with position restraints on Carbon
        Alpha with Low restraints with dt = dt

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param pdb_restr: reference coordinate file for position restraints
        :type pdb_restr: str, default=None

        :param nsteps: number of equilibration steps with BB constraints
        :type nsteps: int, default=100000

        :param dt: integration time step for BB equilibration
        :type dt: float, default=0.002

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if name is None:
            name = self.name
        if pdb_restr is None:
            pdb_restr = self.coor_file

        equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                         "template/martini_v2.x_new-rf.mdp")

        mdp_options.update({'nsteps': int(nsteps),
                            'define': '-DPOSRES',
                            'dt': dt})
        self.run_md_sim(out_folder=out_folder, name="equi_CG_BB_" + name,
                        pdb_restr=pdb_restr, mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

    def prod_CG(self, out_folder, name=None, nsteps=5000000, dt=0.02,
                maxwarn=0,
                monitor_tool=monitor.PROGRESS_BAR,
                **mdp_options):
        """Equilibrate a system a CG system:

        1. equilibration of nsteps_HA with position restraints on Heavy Atoms
        with dt = dt_HA
        2. equilibration of nsteps_CA with position restraints on Carbon Alpha
        with dt = dt
        3. equilibration of nsteps_CA_LOW with position restraints on Carbon
        Alpha with Low restraints with dt = dt

        :param out_folder: path of the output file folder
        :type out_folder: str

        :param name: name of the simulation to run
        :type name: str, default=None

        :param nsteps: number of equilibration steps with BB constraints
        :type nsteps: int, default=100000

        :param dt: integration time step for BB equilibration
        :type dt: float, default=0.002

        :param monitor: option to monitor a simulation, if not none monitor
            should contains two values: ``function`` the function to be ran
            while simulation is running and ``input`` parameters for the
            function
        :type monitor: dict, default=None

        :param mdp_options: Additional mdp parameters to use
        :type mdp_options: dict

        **Object requirement(s):**

            * self.coor_file
            * self.top_file
            * self.nt
            * self.ntmpi
            * self.gpu_id

        **Object field(s) changed:**

            * self.mdp
            * self.tpr
            * self.coor_file
            * self.xtc
        """

        if name is None:
            name = self.name

        equi_template_mdp = os.path.join(GROMACS_MOD_DIRNAME,
                                         "template/martini_v2.x_new-rf.mdp")

        mdp_options.update({'nsteps': int(nsteps), 'define': '', 'dt': dt})
        self.run_md_sim(out_folder=out_folder, name="prod_CG_" + name,
                        mdp_template=equi_template_mdp,
                        mdp_options=mdp_options, maxwarn=maxwarn,
                        monitor_tool=monitor_tool)

    def get_last_output(self):
        """In a case of a simulation restart, outputs edr, log, gro and xtc
        files are called for example as  ``self.sim_name+".partXXXX.edr"``
        where XXXX is the iteration number of restart (eg. first restart:
        XXXX=0002).

        This function actualise the edr, log, coor_file and xtc variable.

        **Object requirement(s):**

            * self.sim_name

        **Object field(s) changed:**

            * self.edr
            * self.log
            * self.coor_file
            * self.xtc

        """

        import glob
        # Get all edr files name :
        edr_file_list = glob.glob(self.sim_name + '.part*.edr')

        index_list = [int(file[-8:-4]) for file in edr_file_list]

        last_index = max(index_list)
        self.edr = '{}.part{:04d}.edr'.format(self.sim_name, last_index)
        self.log = '{}.part{:04d}.log'.format(self.sim_name, last_index)
        self.coor_file = '{}.part{:04d}.gro'.format(self.sim_name, last_index)
        self.xtc = '{}.part{:04d}.xtc'.format(self.sim_name, last_index)
        if os_command.check_file_exist(
                '{}.part{:04d}.xvg'.format(self.sim_name, last_index)):
            self.xvg = '{}.part{:04d}.xvg'.format(self.sim_name, last_index)

    def get_all_output(self):
        """In a case of a simulation restart, outputs edr, log, gro, xvg and
        xtc files are called for example as  ``self.sim_name+".partXXXX.edr"``
        where XXXX is the iteration number of restart (eg. first restart:
        XXXX=0002).

        This function return a dictionnary of all edr, log, coor_file, xvg
        and xtc list.

        **Object requirement(s):**

            * self.sim_name

        :return: return dict containing edr, log, xtc, xvg and coor_file
            file list
        :rtype: dict

        """

        import glob
        # Get all edr files name :
        edr_file_list = glob.glob(self.tpr[:-4] + '*.edr')

        index_list = [file[:-4] for file in edr_file_list]

        output_dict = {'edr': [], 'log': [], 'gro': [], 'xtc': [], 'xvg': []}

        for index in index_list:
            output_dict['edr'].append('{}.edr'.format(index))
            output_dict['log'].append('{}.log'.format(index))
            output_dict['gro'].append('{}.gro'.format(index))
            output_dict['xtc'].append('{}.xtc'.format(index))
            output_dict['xvg'].append('{}.xvg'.format(index))

        return(output_dict)

    def get_simulation_time(self):
        """In a case of a simulation restart simulation, one would like to
        know how much simulation time has already been computed to reach a
        certain amount of time in the restart simulation.
        The command will check the cpt file using ``gmx check -f file.cpt``.

        **Object requirement(s):**

            * self.tpr

        **Object field(s) changed:**

            * None

        :return: return simulation time (ns)
        :rtype: float

        """

        cpt_file = self.tpr[:-4] + ".cpt"
        logger.info("- Get simulation time from : {}".format(cpt_file))

        # Check if output files exist:
        if not os.path.isfile(cpt_file):
            logger.info("Checkpoint file {} coulnd not be found".format(
                cpt_file))
            raise FileNotFoundError()

        cmd_list = [GMX_BIN, "check",
                    "-f", cpt_file]

        cmd_run = os_command.Command(cmd_list)

        cmd_run.display()
        output = cmd_run.run(out_data=True)

        # Search in all line, if it start with "Last frame"
        for line in output['stderr'].splitlines():
            if line.startswith("Last frame"):
                time = float(line.split()[4])
                return time

        logger.error("Last Frame not found in gmx check output")
        raise ValueError()

    ##########################################################
    # ###########   ANALYSIS RELATED FUNCTIONS   #############
    ##########################################################

    def convert_selection_to_index(self, selection_list):
        """ Convert selection list with selection name eg.
        ["System"] to the index number eg. ["0"].

        :param selection_list: List of selection names
        :type selection_list: list

        :return: list of selection numbers
        :rtype: list

        """

        sele_dict = self.get_index_dict()

        sele_index_list = []
        for sele in selection_list:
            if sele in sele_dict:
                sele_index_list.append(str(sele_dict[sele]))
            else:
                try:
                    int(sele)
                except ValueError:
                    logger.warning('Selection {} could not be found'
                                   ' in the index file {}'.format(
                                        sele, self.ndx))
        return(sele_index_list)

    def get_ener(self, selection_list, output_xvg='tmp_edr.xvg',
                 check_file_out=True, keep_ener_file=False):
        """Get energy of a system using ``gmx energy``.

        :param selection_list: List of selection names or number
        :type selection_list: list

        :param output_xvg: output `.xvg` file name
        :type output_xvg: str, default='tmp_edr.xvg'

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, default=True

        :param keep_ener_file: flag to keep or not output `.xvg` file.
        :type keep_ener_file: bool, default=False

        :return: Energy table
        :rtype: pd.DataFrame

       """

        logger.info("- Extract energy")

        # Check if output files exist:
        if check_file_out and os.path.isfile(output_xvg):
            logger.info("get_ener not launched {} already exist".format(
                output_xvg))
        else:
            cmd_convert = os_command.Command([GMX_BIN, "energy",
                                              "-f", self.edr,
                                              "-o", output_xvg])

            cmd_convert.display()
            cmd_convert.run(com_input='\n'.join(selection_list))

        ener_pd = monitor.read_xvg(output_xvg)

        if not keep_ener_file:
            os_command.delete_file(output_xvg)

        return(ener_pd)

    def get_rmsd(self, selection_list=['C-alpha', 'Protein'],
                 output_xvg='tmp_rmsd.xvg',
                 fit="rot+trans", pbc="no", ref_sys=None,
                 keep_ener_file=False):
        """Get RMSD of a system using ``gmx rms``.


        :param selection_list: List of selection names or number
        :type selection_list: list, default=['C-alpha', 'Protein']

        :param output_xvg: output `.xvg` file name
        :type output_xvg: str, default='tmp_rmsd.xvg'

        :param fit: Coordinates Fitting Method
        :type fit: str, default="rot+trans"

        :param pbc: Peridic Boundary condition treatment
        :type pbc: str, default="no"

        :param keep_ener_file: flag to keep or not output `.xvg` file.
        :type keep_ener_file: bool, default=False

        :param ref_sys: reference system for RMSD calculation
        :type ref_sys: GmxSys, default=None

        :return: RMSD table
        :rtype: pd.DataFrame
        """

        logger.info("- Extract RMSD")

        if ref_sys is None:
            ref_sys = self

        # Convert selection name to index:
        sele_index_list = self.convert_selection_to_index(selection_list)
        group_num = len(sele_index_list) - 1

        cmd_list = [GMX_BIN, "rms",
                    "-s", ref_sys.tpr,
                    "-f", self.xtc,
                    "-n", self.ndx,
                    "-o", output_xvg,
                    "-fit", fit,
                    "-ng", str(group_num),
                    "-pbc", pbc]

        cmd_rmsd = os_command.Command(cmd_list)

        cmd_rmsd.display()
        cmd_rmsd.run(com_input='\n'.join(sele_index_list))

        ener_pd = monitor.read_xvg(output_xvg)
        ener_pd.columns = ['time'] + selection_list[1:]
        if not keep_ener_file:
            os_command.delete_file(output_xvg)

        return(ener_pd)

    def get_rmsf(self, selection_list, output_xvg='tmp_rmsf.xvg',
                 fit="no", res="no",
                 keep_ener_file=False):
        """Get RMSF of a system using ``gmx rmsf``.

        :param selection_list: List of selection names or number
        :type selection_list: list

        :param output_xvg: output `.xvg` file name
        :type output_xvg: str, default='tmp_rmsf.xvg'

        :param fit: Flag for fitting before computing RMSF
        :type fit: str, default="no"

        :param res: Residue averaging flag
        :type res: str, default="no"

        :param keep_ener_file: flag to keep or not output `.xvg` file.
        :type keep_ener_file: bool, default=False

        :return: RMSF table
        :rtype: pd.DataFrame
        """

        logger.info("- Extract RMSF")

        # Convert selection name to index:
        sele_index_list = self.convert_selection_to_index(selection_list)

        cmd_list = [GMX_BIN, "rmsf",
                    "-s", self.tpr,
                    "-f", self.xtc,
                    "-n", self.ndx,
                    "-o", output_xvg,
                    "-fit", fit,
                    "-res", res]

        cmd_rmsd = os_command.Command(cmd_list)

        cmd_rmsd.display()
        cmd_rmsd.run(com_input='\n'.join(sele_index_list))

        ener_pd = monitor.read_xvg(output_xvg)
        ener_pd.rename(columns={'(nm)': 'RMSF'},
                       inplace=True)

        if not keep_ener_file:
            os_command.delete_file(output_xvg)

        return(ener_pd)

    def get_dist(self, distance_list, output_xvg='tmp_dist.xvg',
                 keep_ener_file=False):
        """Get distances as a function of time
        for a trajectory using ``gmx distance``.


        :param distance_list: List of atom couple
        :type distance_list: list

        :param output_xvg: output `.xvg` file name
        :type output_xvg: str, default='tmp_dist.xvg'

        :param keep_ener_file: flag to keep or not output `.xvg` file.
        :type keep_ener_file: bool, default=False

        :return: Distance table
        :rtype: pd.DataFrame
        """

        logger.info("- Extract distance")

        # Write index with all distances pairs:
        with open('tmp.ndx', 'w') as file_out:
            for i, dist in enumerate(distance_list):
                file_out.write('[ dist_{} ]\n'.format(i))
                file_out.write('{} {}\n'.format(dist[0], dist[1]))

        cmd_convert = os_command.Command([GMX_BIN, "distance",
                                          "-n", 'tmp.ndx',
                                          "-f", self.xtc,
                                          "-s", self.tpr,
                                          "-oall", output_xvg])

        cmd_convert.display()
        # cmd_convert.run(com_input='\n'.join(selection_list))
        cmd_convert.run(com_input='\n'.join(
            [str(i) for i in range(len(distance_list))]))

        ener_pd = monitor.read_xvg(output_xvg)
        # print(ener_pd)
        # ener_pd = ener_pd.reset_index()
        ener_pd.columns = [ener_pd.columns[0]] +\
            ['dist_{}_{}'.format(dist[0], dist[1]) for dist in distance_list]

        if not keep_ener_file:
            os_command.delete_file(output_xvg)
        os_command.delete_file('tmp.ndx')

        return(ener_pd)

    def get_angle(self, angle_list, output_xvg='tmp_angle.xvg',
                  keep_ener_file=False, improper=False):
        """Get angle of a traj using ``gmx angle``.

        :param angle_list: List of atom triplet
        :type angle_list: list

        :param output_xvg: output `.xvg` file name
        :type output_xvg: str, default='tmp_angle.xvg'

        :param keep_ener_file: flag to keep or not output `.xvg` file.
        :type keep_ener_file: bool, default=False

        :param improper: flag to compute improper angles.
        :type improper: bool, default=False

        :return: Distance table
        :rtype: pd.DataFrame

        """

        # logger.info("- Extract angle")
        if improper:
            angle_type = 'improper'
        elif len(angle_list[0]) == 3:
            angle_type = 'angle'
        elif len(angle_list[0]) == 4:
            angle_type = 'dihedral'
        else:
            print('Error')

        # Write index with all angle triplets/quadruplets:
        label_list = []
        with open('tmp.ndx', 'w') as file_out:
            file_out.write('[ angles ]\n')
            for i, angle in enumerate(angle_list):
                angle_str = ''
                for j in angle:
                    angle_str += '{} '.format(j)
                file_out.write('{}\n'.format(angle_str))
                label_list.append('angle_' + angle_str[:-1].replace(" ", "_"))

        output_distri_xvg = 'tmp_angdist.xvg'

        cmd_convert = os_command.Command([GMX_BIN, "angle",
                                          "-n", 'tmp.ndx',
                                          "-f", self.xtc,
                                          "-ov", output_xvg,
                                          "-od", output_distri_xvg,
                                          "-type", angle_type,
                                          "-all"])

        cmd_convert.display()
        cmd_convert.run()

        ener_pd = monitor.read_xvg(output_xvg)
        # ener_pd = ener_pd.reset_index()
        # print(ener_pd)
        ener_pd.columns = [ener_pd.columns[2], 'Avg'] +\
            [label for label in label_list]
        ener_pd = ener_pd.drop(['Avg'], axis=1)

        if not keep_ener_file:
            os_command.delete_file(output_xvg)
        os_command.delete_file(output_distri_xvg)
        os_command.delete_file('tmp.ndx')

        return(ener_pd)


if __name__ == "__main__":

    import doctest
    import shutil

    TEST_DIR = 'gromacs_py_test_out'
    TEST_OUT = os.path.join(TEST_DIR, 'gmx5')

    def getfixture(*args):
        return TEST_OUT

    print("-Test gmx5 module:")
    print("gmx5:    \t", doctest.testmod())

    # Erase all test files
    shutil.rmtree(TEST_DIR, ignore_errors=True)
