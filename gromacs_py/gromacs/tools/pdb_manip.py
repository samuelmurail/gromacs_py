#!/usr/bin/env python3
# coding: utf-8

"""
#####################################
#########    PDB IN/OUT    ##########
#####################################
"""

__author__ = "Samuel Murail"

import os
import sys
import time
import numpy as np
from numpy.linalg import norm
from numpy import sin, cos
from scipy.spatial import distance_matrix

# In case pdb_manip is launched as main, relative import will failed
try:
    from . import os_command
except ImportError:
    print("Relative import from . fails, use absolute import instead")
    import os_command

# Test folder path
PDB_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(PDB_LIB_DIR, "../../test/input/"))

# Global variables:
AA_DICT = {'GLY': 'G',
           'HIS': 'H',
           'HSE': 'H',
           'HSD': 'H',
           'HSP': 'H',
           'ARG': 'R',
           'LYS': 'K',
           'ASP': 'D',
           'GLU': 'E',
           'SER': 'S',
           'THR': 'T',
           'ASN': 'N',
           'GLN': 'Q',
           'CYS': 'C',
           'SEC': 'U',
           'PRO': 'P',
           'ALA': 'A',
           'ILE': 'I',
           'PHE': 'F',
           'TYR': 'Y',
           'TRP': 'W',
           'VAL': 'V',
           'LEU': 'L',
           'MET': 'M'}

AA_1_TO_3_DICT = {'G': 'GLY',
                  'H': 'HIS',
                  'R': 'ARG',
                  'K': 'LYS',
                  'D': 'ASP',
                  'E': 'GLU',
                  'S': 'SER',
                  'T': 'THR',
                  'N': 'ASN',
                  'Q': 'GLN',
                  'C': 'CYS',
                  'U': 'SEC',
                  'P': 'PRO',
                  'A': 'ALA',
                  'I': 'ILE',
                  'F': 'PHE',
                  'Y': 'TYR',
                  'W': 'TRP',
                  'V': 'VAL',
                  'L': 'LEU',
                  'M': 'MET',
                  'X': 'ACE'}

# Atom names for each residues
BACK_ATOM = ['N', 'CA', 'C', 'O']

AA_ATOM_DICT = {'X': ['CH3', 'O', 'C'],  # X:ACE
                'G': BACK_ATOM,
                'A': BACK_ATOM + ['CB'],
                'S': BACK_ATOM + ['CB', 'OG'],
                'C': BACK_ATOM + ['CB', 'SG'],
                'T': BACK_ATOM + ['CB', 'OG1', 'CG2'],
                'V': BACK_ATOM + ['CB', 'CG1', 'CG2'],
                'I': BACK_ATOM + ['CB', 'CG1', 'CG2', 'CD'],
                'L': BACK_ATOM + ['CB', 'CG', 'CD1', 'CD2'],
                'N': BACK_ATOM + ['CB', 'CG', 'ND2', 'OD1'],
                'D': BACK_ATOM + ['CB', 'CG', 'OD1', 'OD2'],
                'M': BACK_ATOM + ['CB', 'CG', 'SD', 'CE'],
                'Q': BACK_ATOM + ['CB', 'CG', 'CD', 'NE2', 'OE1'],
                'E': BACK_ATOM + ['CB', 'CG', 'CD', 'OE1', 'OE2'],
                'K': BACK_ATOM + ['CB', 'CG', 'CD', 'CE', 'NZ'],
                'R': BACK_ATOM + ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                'F': BACK_ATOM + ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2'],
                'Y': BACK_ATOM + ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2', 'OH'],
                'H': BACK_ATOM + ['CB', 'CG', 'ND1', 'CE1', 'CD2', 'NE2'],
                'W': BACK_ATOM + ['CB', 'CG', 'CD1', 'NE1', 'CD2', 'CE2', 'CE3', 'CZ3', 'CH2', 'CZ2'],
                'P': BACK_ATOM + ['CB', 'CG', 'CD']}

# Bond definition:
# Note that order is important
BACK_BOND = [['-C', 'N'], ['N', 'CA'], ['CA', 'C'], ['C', 'O']]
# X is for ACE special case
AA_BOND_DICT = {}
AA_BOND_DICT['X'] = [['CH3', 'O'], ['O', 'C']]  # Need to use a trick with unphysical bond
AA_BOND_DICT['G'] = BACK_BOND
AA_BOND_DICT['A'] = BACK_BOND + [['CA', 'CB']]
AA_BOND_DICT['S'] = BACK_BOND + [['CA', 'CB'], ['CB', 'OG']]
AA_BOND_DICT['C'] = BACK_BOND + [['CA', 'CB'], ['CB', 'SG']]
AA_BOND_DICT['T'] = BACK_BOND + [['CA', 'CB'], ['CB', 'OG1'], ['CB', 'CG2']]
AA_BOND_DICT['V'] = BACK_BOND + [['CA', 'CB'], ['CB', 'CG1'], ['CB', 'CG2']]
AA_BOND_DICT['I'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG1'], ['CG1', 'CD'], ['CB', 'CG2']]
AA_BOND_DICT['L'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CG', 'CD2']]
AA_BOND_DICT['N'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'ND2'], ['CG', 'OD1']]
AA_BOND_DICT['D'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'OD1'], ['CG', 'OD2']]
AA_BOND_DICT['M'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'SD'], ['SD', 'CE']]
AA_BOND_DICT['Q'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'NE2'],
     ['CD', 'OE1']]
AA_BOND_DICT['E'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'OE1'],
     ['CD', 'OE2']]
AA_BOND_DICT['K'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'CE'],
     ['CE', 'NZ']]
AA_BOND_DICT['R'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'NE'],
     ['NE', 'CZ'], ['CZ', 'NH1'], ['CZ', 'NH2']]
AA_BOND_DICT['F'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'CE1'],
     ['CE1', 'CZ'], ['CG', 'CD2'], ['CD2', 'CE2']]
AA_BOND_DICT['Y'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'CE1'],
     ['CE1', 'CZ'], ['CZ', 'OH'], ['CG', 'CD2'], ['CD2', 'CE2']]
AA_BOND_DICT['H'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'ND1'], ['ND1', 'CE1'],
     ['CG', 'CD2'], ['CD2', 'NE2']]
AA_BOND_DICT['W'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'NE1'],
     ['CG', 'CD2'], ['CD2', 'CE2'], ['CD2', 'CE3'], ['CE3', 'CZ3'],
     ['CZ3', 'CH2'], ['CH2', 'CZ2']]
AA_BOND_DICT['P'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD']]

# Distance, angle and dihedral angles parameters
BACK_DIST = [['N', 'CA', 1.46],
             ['CA', 'C', 1.52],
             ['C', 'O', 1.23],
             ['C', 'N', 1.29]]

BACK_ANGLE = [['N', 'CA', 'C', 110.9],
              ['CA', 'C', 'O', 122.0],
              ['CA', 'C', 'N', 110.9],
              ['CA', 'N', 'C', 121.3],
              ['N', 'C', 'O', 119.0]]  # Only for ACE-connexion

BACK_DIHE = [['N', 'CA', 'C', 'O', 0],
             ['N', 'CA', 'C', 'N', 180.0],
             ['CA', 'N', 'C', 'CA', 180.0],
             ['C', 'CA', 'N', 'C', -180.0],
             ['N', 'C', 'O', 'CH3', 180.0],  # Only for ACE-connexion
             ['CA', 'N', 'C', 'O', 0.0]]  # Only for ACE-connexion

DIST_DICT = {}
ANGLE_DICT = {}
DIHE_DICT = {}

# ACE X
DIST_DICT['X'] = [['CH3', 'O', 2.40], ['C', 'O', 1.23]]
ANGLE_DICT['X'] = [['CH3', 'O', 'C', 32.18]]
DIHE_DICT['X'] = BACK_DIHE


# Glycine
DIST_DICT['G'] = BACK_DIST
ANGLE_DICT['G'] = BACK_ANGLE
DIHE_DICT['G'] = BACK_DIHE

# Alanine
DIST_DICT['A'] = BACK_DIST + [['CA', 'CB', 1.52]]
ANGLE_DICT['A'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7]]
DIHE_DICT['A'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3]]

# Serine
DIST_DICT['S'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'OG', 1.42]]
ANGLE_DICT['S'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'OG', 110.8]]
DIHE_DICT['S'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'OG', 69.4]]
# Cysteine
DIST_DICT['C'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'SG', 1.81]]
ANGLE_DICT['C'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'SG', 110.8]]
DIHE_DICT['C'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'SG', -173.8]]
# Threonine
DIST_DICT['T'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'OG1', 1.42],
                              ['CB', 'CG2', 1.54]]
ANGLE_DICT['T'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'OG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3]]
DIHE_DICT['T'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'OG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6]]

# Valine
DIST_DICT['V'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG1', 1.54],
                              ['CB', 'CG2', 1.54]]
ANGLE_DICT['V'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3]]
DIHE_DICT['V'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6]]

# Isoleucine
DIST_DICT['I'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG1', 1.54],
                              ['CB', 'CG2', 1.54],
                              ['CG1', 'CD', 1.54]]
ANGLE_DICT['I'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3],
                                ['CB', 'CG1', 'CD', 116.3]]
DIHE_DICT['I'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6],
                              ['CA', 'CB', 'CG1', 'CD', 179.6]]

# Isoleucine
DIST_DICT['L'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.54],
                              ['CG', 'CD2', 1.54]]
ANGLE_DICT['L'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 110.6],
                                ['CB', 'CG', 'CD2', 116.3]]
DIHE_DICT['L'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -57.8],
                              ['CA', 'CB', 'CG', 'CD1', -61.5],
                              ['CA', 'CB', 'CG', 'CD2', 179.6]]

# Asparagine
DIST_DICT['N'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'ND2', 1.29],
                              ['CG', 'OD1', 1.23]]
ANGLE_DICT['N'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'ND2', 118.9],
                                ['CB', 'CG', 'OD1', 122.2]]
DIHE_DICT['N'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -57.8],
                              ['CA', 'CB', 'CG', 'ND2', -78.2],
                              ['CA', 'CB', 'CG', 'OD1', 100.6]]

# Aspartic Acid
DIST_DICT['D'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'OD1', 1.23],
                              ['CG', 'OD2', 1.23]]
ANGLE_DICT['D'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'OD1', 118.9],
                                ['CB', 'CG', 'OD2', 122.2]]
DIHE_DICT['D'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -177.0],
                              ['CA', 'CB', 'CG', 'OD1', 37.0],
                              ['CA', 'CB', 'CG', 'OD2', -140.7]]

# Methionine
DIST_DICT['M'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'SD', 1.80],
                              ['SD', 'CE', 1.80]]
ANGLE_DICT['M'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'SD', 118.9],
                                ['CG', 'SD', 'CE', 98.5]]
DIHE_DICT['M'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -68.5],
                              ['CA', 'CB', 'CG', 'SD', -165.1],
                              ['CB', 'CG', 'SD', 'CE', -140.7]]

# Glutamine
DIST_DICT['Q'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'NE2', 1.31],
                              ['CD', 'OE1', 1.22]]
ANGLE_DICT['Q'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'NE2', 121.1],
                                ['CG', 'CD', 'OE1', 120.0]]
DIHE_DICT['Q'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -165.1],
                              ['CB', 'CG', 'CD', 'NE2', -0.9],
                              ['CB', 'CG', 'CD', 'OE1', 177.9]]

# Glutamic Acid
DIST_DICT['E'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'OE1', 1.22],
                              ['CD', 'OE2', 1.22]]
ANGLE_DICT['E'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'OE1', 121.1],
                                ['CG', 'CD', 'OE2', 120.0]]
DIHE_DICT['E'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'OE1', -0.9],
                              ['CB', 'CG', 'CD', 'OE2', 177.9]]

# Lysine
DIST_DICT['K'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'CE', 1.54],
                              ['CE', 'NZ', 1.3]]
ANGLE_DICT['K'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'CE', 118.9],
                                ['CD', 'CE', 'NZ', 109.4]]
DIHE_DICT['K'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'CE', 177.1],
                              ['CG', 'CD', 'CE', 'NZ', 177.9]]

# Arginine
DIST_DICT['R'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'NE', 1.54],
                              ['NE', 'CZ', 1.3],
                              ['NH1', 'CZ', 1.3],
                              ['NH2', 'CZ', 1.3]]
ANGLE_DICT['R'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'NE', 118.9],
                                ['CD', 'NE', 'CZ', 125.3],
                                ['NE', 'CZ', 'NH1', 123.6],
                                ['NE', 'CZ', 'NH2', 123.6]]
DIHE_DICT['R'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'NE', 177.1],
                              ['CG', 'CD', 'NE', 'CZ', 177.9],
                              ['CD', 'NE', 'CZ', 'NH1', 0.3],
                              ['CD', 'NE', 'CZ', 'NH2', -179.6]]

# Phenylalanine
DIST_DICT['F'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'CE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CE1', 'CZ', 1.4]]
ANGLE_DICT['F'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 120.8],
                                ['CB', 'CG', 'CD2', 120.8],
                                ['CG', 'CD1', 'CE1', 120.4],
                                ['CG', 'CD2', 'CE2', 120.4],
                                ['CD1', 'CE1', 'CZ', 120.4]]
DIHE_DICT['F'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -65.4],
                              ['CA', 'CB', 'CG', 'CD1', -78.1],
                              ['CA', 'CB', 'CG', 'CD2', 101.4],
                              ['CB', 'CG', 'CD1', 'CE1', 179.1],
                              ['CB', 'CG', 'CD2', 'CE2', 179.1],
                              ['CG', 'CD1', 'CE1', 'CZ', 0.1]]

# Tyrosine
DIST_DICT['Y'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'CE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CE1', 'CZ', 1.4],
                              ['CZ', 'OH', 1.22]]
ANGLE_DICT['Y'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 120.8],
                                ['CB', 'CG', 'CD2', 120.8],
                                ['CG', 'CD1', 'CE1', 120.4],
                                ['CG', 'CD2', 'CE2', 120.4],
                                ['CD1', 'CE1', 'CZ', 120.4],
                                ['CE1', 'CZ', 'OH', 120.8]]
DIHE_DICT['Y'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -65.4],
                              ['CA', 'CB', 'CG', 'CD1', -78.1],
                              ['CA', 'CB', 'CG', 'CD2', 101.4],
                              ['CB', 'CG', 'CD1', 'CE1', 179.1],
                              ['CB', 'CG', 'CD2', 'CE2', 179.1],
                              ['CG', 'CD1', 'CE1', 'CZ', 0.1],
                              ['CD1', 'CE1', 'CZ', 'OH', 179.1]]

# Histidine
DIST_DICT['H'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'ND1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['ND1', 'CE1', 1.4],
                              ['CD2', 'NE2', 1.4]]
ANGLE_DICT['H'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'ND1', 131.5],
                                ['CB', 'CG', 'CD2', 117.8],
                                ['CG', 'ND1', 'CE1', 105.4],
                                ['CG', 'CD2', 'NE2', 105.7]]
DIHE_DICT['H'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -109.4],
                              ['CA', 'CB', 'CG', 'ND1', 107.6],
                              ['CA', 'CB', 'CG', 'CD2', -75.6],
                              ['CB', 'CG', 'ND1', 'CE1', 177.5],
                              ['CB', 'CG', 'CD2', 'NE2', 171.7]]

# Tryptophan
DIST_DICT['W'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'NE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CD2', 'CE3', 1.4],
                              ['CE3', 'CZ3', 1.4],
                              ['CZ3', 'CH2', 1.4],
                              ['CH2', 'CZ2', 1.4]]
ANGLE_DICT['W'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 131.5],
                                ['CB', 'CG', 'CD2', 117.8],
                                ['CG', 'CD1', 'NE1', 105.4],
                                ['CG', 'CD2', 'CE2', 105.7],
                                ['CG', 'CD2', 'CE3', 131.5],
                                ['CD2', 'CE3', 'CZ3', 120.4],
                                ['CE3', 'CZ3', 'CH2', 120.4],
                                ['CZ3', 'CH2', 'CZ2', 120.4]]
DIHE_DICT['W'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -109.4],
                              ['CA', 'CB', 'CG', 'CD1', 107.6],
                              ['CA', 'CB', 'CG', 'CD2', -75.6],
                              ['CB', 'CG', 'CD1', 'NE1', 177.5],
                              ['CB', 'CG', 'CD2', 'CE2', 180],
                              ['CB', 'CG', 'CD2', 'CE3', 0.0],
                              ['CG', 'CD2', 'CE3', 'CZ3', 180.0],
                              ['CD2', 'CE3', 'CZ3', 'CH2', 0.0],
                              ['CE3', 'CZ3', 'CH2', 'CZ2', 0]]

# Proline
DIST_DICT['P'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54]]
ANGLE_DICT['P'] = BACK_ANGLE + [['CB', 'CA', 'N', 101.9],
                                ['CA', 'CB', 'CG', 103.7],
                                ['CB', 'CG', 'CD', 103.3]]

DIHE_DICT['P'] = [['N', 'CA', 'C', 'O', 0],
                  ['N', 'CA', 'C', 'N', 180.0],
                  ['CA', 'N', 'C', 'CA', 180.0],
                  ['C', 'CA', 'N', 'C', -70.0],
                  ['N', 'C', 'O', 'CH3', 180.0],  # Only for ACE-connexion
                  ['CA', 'N', 'C', 'O', 0.0],
                  ['CB', 'CA', 'N', 'C', 168.6],
                  ['N', 'CA', 'CB', 'CG', 29.6],
                  ['CA', 'CB', 'CG', 'CD', -37.5]]


class Coor:
    """ Topologie base on coordinates like pdb or gro.

    The coor object containt a dictionnary of atoms indexed
    on the atom num and the crystal packing info.


    :param atom_dict: dictionnary of atom
    :type atom_dict: dict

    :param crystal_pack: crystal packing
    :type crystal_pack: str

    **Atom dictionnary parameters**

    :param field: pdb field
    :type field: str

    :param num: atom number
    :type num: int

    :param name: atom name
    :type name: str

    :param alter_loc: atom number
    :type alter_loc: str

    :param res_name: residue name (3 letters)
    :type res_name: str

    :param chain: chain ID
    :type chain: str

    :param res_num: residue number (based on pdb file)
    :type res_num: int

    :param uniq_resid: unique residue number
    :type uniq_resid: int

    :param insert_res: atom number
    :type insert_res: str

    :param xyz: coordinate
    :type x: numpy array

    :param occ: occupation
    :type occ: float

    :param beta: beta flactor
    :type beta: float


    .. note::
        The atom num index in the dictionnary, is not the same as the
        ``atom_num`` field of the dictionnary.

    .. note::
        Files necessary for testing : ../test/input/1y0m.pdb, ../test/input/1rxz.pdb
        and ../test/input/4n1m.pdb.
        To do the unitary test, execute pdb_mani.py (-v for verbose mode)

    .. todo::
        Add an atom class ?

    """

    def __init__(self):
        self.atom_dict = dict()
        self.crystal_pack = None

    def read_pdb(self, pdb_in, pqr_format=False):
        """Read a pdb file and return atom informations as a dictionnary indexed on the atom num.
        The fonction can also read pqr files if specified with ``pqr_format = True``,
        it will only change the column format of beta and occ factors.

        :param pdb_in: path of the pdb file to read
        :type pdb_in: str

        :param pqr_format: Flag for .pqr file format reading.
        :type pqr_format: bool, default=False

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found

        """

        atom_index = 0
        uniq_resid = -1
        old_res_num = -1

        with open(pdb_in) as pdbfile:
            for line in pdbfile:
                if line[:6] == "CRYST1":
                    self.crystal_pack = line
                if line[:4] == 'ATOM' or line[:6] == "HETATM":

                    field = line[:6].strip()
                    atom_num = int(line[6:11])
                    atom_name = line[12:16].strip()
                    
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26])
                    insert_res = line[26:27]
                    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    if pqr_format:
                        alter_loc = ""
                        res_name = line[16:20].strip()
                        occ, beta = line[54:62].strip(), line[62:70].strip()
                    else:
                        alter_loc = line[16:17]
                        res_name = line[17:20].strip()
                        occ, beta = line[54:60].strip(), line[60:66].strip()

                    if occ == "":
                        occ = 0.0
                    else:
                        occ = float(occ)

                    if beta == "":
                        beta = 0.0
                    else:
                        beta = float(beta)

                    if res_num != old_res_num:
                        uniq_resid += 1
                        old_res_num = res_num

                    atom = {"field": field,
                            "num": atom_num,
                            "name": atom_name,
                            "alter_loc": alter_loc,
                            "res_name": res_name,
                            "chain": chain,
                            "res_num": res_num,
                            "uniq_resid": uniq_resid,
                            "insert_res": insert_res,
                            "xyz": xyz,
                            "occ": occ,
                            "beta": beta}

                    self.atom_dict[atom_index] = atom
                    atom_index += 1
        print("Succeed to read file", os.path.relpath(pdb_in), ", ", atom_index, "atoms found")

    def write_pdb(self, pdb_out, check_file_out=True):
        """Write a pdb file.

        :param pdb_out: path of the pdb file to write
        :type pdb_out: str

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.write_pdb(os.path.join(TEST_OUT, 'tmp.pdb')) #doctest: +ELLIPSIS
        Succeed to save file .../tmp.pdb

        """

        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            print("PDB file {} already exist, file not saved".format(pdb_out))
            return

        filout = open(pdb_out, 'w')
        if self.crystal_pack is not None:
            filout.write(self.crystal_pack)
        for atom_num, atom in sorted(self.atom_dict.items()):
            # print(pdb_dict[atom_num]["name"])
            filout.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                atom["field"],
                atom["num"],
                atom["name"],
                atom["alter_loc"],
                atom["res_name"],
                atom["chain"],
                atom["res_num"],
                atom["insert_res"],
                atom["xyz"][0],
                atom["xyz"][1],
                atom["xyz"][2],
                atom["occ"],
                atom["beta"]))
        filout.write("TER\n")
        filout.close()

        print("Succeed to save file", os.path.relpath(pdb_out))
        return

    def get_aa_seq(self):
        """Get the amino acid sequence from a coor object.

        :return: dictionnary of sequence indexed by the chain ID
        :rtype: dict

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        .. warning::
            If atom chains are not arranged sequentialy (A,A,A,B,B,A,A,A ...),
            the first atom seq will be overwritten by the last one.

        """

        # Get CA atoms
        CA_index_list = self.get_index_selection({"name": ["CA"]})

        seq = ""
        seq_dict = {}
        chain_first = self.atom_dict[CA_index_list[0]]['chain']

        for index in sorted(CA_index_list):
            loop_atom = self.atom_dict[index]

            if loop_atom['chain'] != chain_first:
                seq_dict[chain_first] = seq
                seq = ""
                chain_first = loop_atom['chain']

            seq = seq + AA_DICT[loop_atom['res_name']]

        seq_dict[chain_first] = seq

        return seq_dict

    def get_aa_num(self):
        """Get the amino acid number of a coor object.

        :return: Number of residues
        :rtype: int

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61

        .. note::
            Only count Ca atoms, this may not be the best choice ?

        """

        CA_index_list = self.get_index_selection({"name": ["CA"]})

        return len(CA_index_list)

    def change_pdb_field(self, change_dict):
        """Change all atom field of a coor object,
        the change is based on the change_dict dictionnary.

        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.change_pdb_field(change_dict = {"chain" : "B"}) #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'B': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        """

        for atom_num, atom in self.atom_dict.items():
            for change, val in change_dict.items():
                atom[change] = val

        return self

    def change_index_pdb_field(self, index_list, change_dict):
        """Change all atom field of a part of coor object defined by ``index``,
        the change is based on the change_dict dictionnary.

        :param index_list: list of atom index to change
        :type index_list: list

        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> res_826_852 = prot_coor.get_index_selection({'res_num' : range(826,852)})
        >>> prot_coor.change_index_pdb_field(index_list = res_826_852, change_dict = {"chain" : "B"}) #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_seq = prot_coor.get_aa_seq()
        >>> prot_seq == {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQD', 'B': 'GGWWRGDYGGKKQLWFPSNYVEEMIN'}
        True

        """

        for atom_num in index_list:
            for change, val in change_dict.items():
                self.atom_dict[atom_num][change] = val

        return self

    def select_part_dict(self, selec_dict):
        """Select atom of a coor object defined,
        the selection is based on the change_dict dictionnary.
        Return a new coor object.

        :param selec_dict: change ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: a new coor object
        :rtype: coor

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61
        >>> prot_20_coor = prot_coor.select_part_dict(selec_dict = {'res_num' : list(range(791,800))})
        >>> prot_20_coor.get_aa_seq()
        {'A': 'TFKSAVKAL'}
        >>> prot_20_coor.get_aa_num()
        9
        >>> prot_N_atom = prot_coor.select_part_dict(selec_dict = {'name' : ['ZN']})
        >>> # WARNING using selec_dict = {'name' : 'ZN'} will give you 61 residues !!
        >>> print(len(prot_N_atom.atom_dict))
        0

        """

        coor_out = Coor()

        for atom_num, atom in self.atom_dict.items():
            selected = True
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                coor_out.atom_dict[atom_num] = atom

        return coor_out

    def get_index_selection(self, selec_dict):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of index of selected atoms.

        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: list of atom index
        :rtype: list of int

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_index_selection({'res_num' : [826,827]})
        [297, 298, 299, 300, 301, 302, 303, 304]

        """

        index_list = []

        for atom_num, atom in self.atom_dict.items():
            selected = True
            # print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("selection=",selection)
                # print("atom:",atom)
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                # print(atom)
                index_list.append(atom_num)

        return index_list

    def get_attribute_selection(self, selec_dict={}, attribute='uniq_resid', index_list=None):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of unique attribute of the selected atoms.

        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: list of atom index
        :rtype: list of int

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_attribute_selection({'res_num' : [826,827]}, attribute='uniq_resid')
        [35, 36]

        """

        attr_list = []

        if index_list is None:
            index_list = self.atom_dict.keys()

        # for atom_num, atom in sorted(self.atom_dict.items()):
        for atom_num in index_list:
            atom = self.atom_dict[atom_num]
            selected = True
            # print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("selection=",selection)
                # print("atom:",atom)
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                # print(atom)
                attr_list.append(atom[attribute])

        return list(set(attr_list))

    def del_atom_index(self, index_list):
        """Delete atoms of a coor object defined by their ``index``.

        :param index_list: list of atom index to delete
        :type index_list: list

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> res_810_852 = prot_coor.get_index_selection({'res_num' : range(810,852)})
        >>> prot_coor.del_atom_index(index_list = res_810_852) #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDE'}

        """

        for index in index_list:
            # print(index, self.atom_dict[index])
            del self.atom_dict[index]

        return self

    def correct_chain(self, Ca_cutoff=4.5):
        """Correct the chain ID's of a coor object, by checking consecutive Calphas atoms distance.
        If the distance is higher than ``Ca_cutoff``, the former atoms are considered
        as in a different chain.

        :param Ca_cutoff: cutoff for distances between Calphas atoms (X)
        :type Ca_cutoff: float, default=4.5

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> res_810 = prot_coor.get_index_selection({'res_num' : [810]})
        >>> prot_coor = prot_coor.del_atom_index(index_list = res_810)
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDETFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.correct_chain() #doctest: +ELLIPSIS
        Chain: A  Residue: 0 to 18
        Chain: B  Residue: 20 to 60
        <...Coor object at ...>
        >>> # As a residue is missing, Calphas after residue 18 is no more consecutive


        .. note::
            This is specially usefull for pdb2gmx which cut the protein
            chains based on the chain ID's.

        """

        Ca_atom = self.select_part_dict({"name": ["CA"]})
        first_flag = True
        chain_res_list = []
        res_list = []

        # Identify Chain uniq_resid
        # Need to use sorted to be sure to check consecutive residues (atoms)
        for key, atom in sorted(Ca_atom.atom_dict.items()):
            if first_flag:
                first_flag = False
            else:
                distance = Coor.atom_dist(atom, old_atom)
                # print(distance)
                if distance < Ca_cutoff:
                    old_atom = atom
                else:
                    # print("New chain")
                    chain_res_list.append(res_list)
                    res_list = []
            res_list.append(atom['uniq_resid'])
            old_atom = atom
        chain_res_list.append(res_list)

        # Change chain ID :

        # print(Ca_atom.atomChabge_dict)

        for i, chain_res in enumerate(chain_res_list):
            print("Chain:", chr(65 + i), " Residue:", chain_res[0], "to", chain_res[-1])
            chain_index = self.get_index_selection({'uniq_resid': chain_res})
            # print(chain_index)
            self.change_index_pdb_field(chain_index, {"chain": chr(65 + i)})
            # pdb_dict_out.update(chain_dict)

        # print(pdb_dict_out)

        return self

    def correct_his_name(self):
        """ Get his protonation state from pdb2pqr and replace HIS resname
        with HSE, HSD, HSP resname.
        To do after pdb2pqr, in order that protonation is recognize by pdb2gmx.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_PATH, '4n1m.pdb'), os.path.join(TEST_OUT, '4n1m.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/4n1m.pdb ,  2530 atoms found
        Succeed to save file .../tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain .../tmp_pdb2pqr.pdb .../4n1m.pqr
        0
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_OUT, '4n1m.pqr'), pqr_format = True) #doctest: +ELLIPSIS
        Succeed to read file .../4n1m.pqr ,  2548 atoms found
        >>> HSD_index = prot_coor.get_index_selection({'res_name' : ['HSD'], 'name':['CA']})
        >>> print(len(HSD_index))
        5
        >>> HSE_index = prot_coor.get_index_selection({'res_name' : ['HSE'], 'name':['CA']})
        >>> print(len(HSE_index))
        0
        >>> HSP_index = prot_coor.get_index_selection({'res_name' : ['HSP'], 'name':['CA']})
        >>> print(len(HSP_index))
        0
        >>> prot_coor.correct_his_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> HIS_index = prot_coor.get_index_selection({'res_name' : ['HIS'], 'name':['CA']})
        >>> print(len(HIS_index))
        0

        .. note::
            This function seems useless. Since last version of pdb2pqr residue name seems correct.
        """

        # FIND HISTIDINE res

        # HSD:
        hsd_uniq_res = self.get_attribute_selection({"res_name": ["HIS"],
                                                     "name": ["HD1"]},
                                                    attribute='uniq_resid')
        # HSE:
        hse_uniq_res = self.get_attribute_selection({"res_name": ["HIS"],
                                                     "name": ["HE2"]},
                                                    attribute='uniq_resid')
        # HSP: find res in common with both hsd and hse
        hsp_uniq_res = [res for res in hsd_uniq_res if res in hse_uniq_res]
        # remove HSP res from HSE HSD list
        if hsp_uniq_res:
            for res in hsp_uniq_res:
                hsd_uniq_res.remove(res)
                hse_uniq_res.remove(res)

        # Replace HIS resname :
        all_his_uniq_res = hsd_uniq_res + hse_uniq_res + hsp_uniq_res

        for atom_num, atom in self.atom_dict.items():
            if atom["uniq_resid"] in all_his_uniq_res:
                if atom["uniq_resid"] in hsd_uniq_res:
                    atom["res_name"] = "HSD"
                elif atom["uniq_resid"] in hse_uniq_res:
                    atom["res_name"] = "HSE"
                else:
                    atom["res_name"] = "HSP"

        return self

    def add_zinc_finger(self, ZN_pdb, cutoff=3.2):
        """ Change protonation state of cysteins and histidine coordinating Zinc atoms.
        To do after `correct_his_name` and `correct_cys_name`, in order that protonation is recognize by pdb2gmx.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>>
        >>> # Read the pdb 1jd4 and keep only chain A
        >>> input_pdb = Coor()
        >>> input_pdb.read_pdb(os.path.join(TEST_PATH, '1jd4.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1jd4.pdb ,  1586 atoms found
        >>> chain_A = input_pdb.select_part_dict(selec_dict = {'chain' : ['A']})
        >>> chain_A.write_pdb(os.path.join(TEST_OUT, '1jd4_A.pdb')) #doctest: +ELLIPSIS
        Succeed to save file .../1jd4_A.pdb
        >>>
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_OUT, '1jd4_A.pdb'), os.path.join(TEST_OUT, '1jd4.pqr')) #doctest: +ELLIPSIS
        Succeed to read file .../1jd4_A.pdb ,  793 atoms found
        Succeed to save file .../tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain .../tmp_pdb2pqr.pdb .../1jd4.pqr
        0
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_OUT, '1jd4.pqr'), pqr_format = True)
        Succeed to read file .../1jd4.pqr ,  1548 atoms found
        >>> prot_coor.correct_cys_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.correct_his_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.correct_chain() #doctest: +ELLIPSIS
        Chain: A  Residue: 0 to 95
        <...Coor object at 0x...
        >>> ZN_index = prot_coor.get_index_selection({'name':['ZN']})
        >>> print(len(ZN_index))
        0
        >>> prot_coor.add_zinc_finger(os.path.join(TEST_OUT, '1jd4_A.pdb')) #doctest: +ELLIPSIS
        Succeed to read file .../1jd4_A.pdb ,  793 atoms found
        Presence of 1 Zinc detected
        change cystein residue(s) : [48, 51, 75]
        change histidine residue(s) : [68]
        True
        >>> ZN_index = prot_coor.get_index_selection({'name':['ZN']})
        >>> print(len(ZN_index))
        1

        .. note::
            This function seems useless. Since last version of pdb2pqr residue name seems correct.
        """

        # Check the number of ZN atoms:
        coor_pre_pqr = Coor()
        coor_pre_pqr.read_pdb(ZN_pdb)
        Zinc_sel = coor_pre_pqr.select_part_dict(selec_dict={'name': ['ZN']})
        Zinc_num = len(Zinc_sel.atom_dict)

        if Zinc_num == 0:
            return False
        else:
            print("Presence of {} Zinc detected".format(Zinc_num))

        # Add the Zinc atoms:
        for key, val in Zinc_sel.atom_dict.items():
            atom = val
            atom['chain'] = 'Z'
            self.atom_dict[len(self.atom_dict)] = val

        # Check cystein and histidine atoms close to ZN:
        close_atom = self.dist_under_index(Zinc_sel, cutoff=cutoff)
        cys_uniq_res_list = []
        his_uniq_res_list = []
        for atom in close_atom:
            local_atom = self.atom_dict[atom]
            if local_atom['res_name'] in ['CYS']:
                cys_uniq_res_list.append(local_atom['uniq_resid'])
            if local_atom['res_name'] in ['HIS', 'HSD', 'HSE', 'HSP']:
                his_uniq_res_list.append(local_atom['uniq_resid'])

        cys_uniq_res_list = list(set(cys_uniq_res_list))
        his_uniq_res_list = list(set(his_uniq_res_list))

        # Change CYS to CYN:
        print("change cystein residue(s) : {}".format(sorted(cys_uniq_res_list)))
        to_change_cys = self.get_index_selection({'uniq_resid': cys_uniq_res_list})
        self.change_index_pdb_field(to_change_cys, {'res_name': 'CYN'})

        # Change Histidine to HSD or HSE
        # the non protonated nitrogen have to be the closest to the ZN atom:
        print("change histidine residue(s) : {}".format(sorted(his_uniq_res_list)))
        for his_uniq_res in his_uniq_res_list:
            # NE2 ND1
            epsilon_his_index = self.get_index_selection({'uniq_resid': [his_uniq_res],
                                                          'name': ['NE2']})[0]
            delta_his_index = self.get_index_selection({'uniq_resid': [his_uniq_res],
                                                        'name': ['ND1']})[0]
            for key, val in Zinc_sel.atom_dict.items():

                epsilon_dist = Coor.atom_dist(val, self.atom_dict[epsilon_his_index])
                delta_dist = Coor.atom_dist(val, self.atom_dict[delta_his_index])

                if (epsilon_dist < cutoff or delta_dist < cutoff):
                    to_change_his = self.get_index_selection({'uniq_resid': [his_uniq_res]})
                    if (epsilon_dist < delta_dist):
                        self.change_index_pdb_field(to_change_his, {'res_name': 'HSD'})
                    else:
                        self.change_index_pdb_field(to_change_his, {'res_name': 'HSE'})

        return True

    def correct_cys_name(self):
        """ Correct the CYS resname from pdb2pqr

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_PATH, '1dpx.pdb'), os.path.join(TEST_OUT, '1dpx.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1dpx.pdb ,  1192 atoms found
        Succeed to save file .../tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain .../tmp_pdb2pqr.pdb .../1dpx.pqr
        0
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_OUT, '1dpx.pqr'), pqr_format = True)
        Succeed to read file .../1dpx.pqr ,  1960 atoms found
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['DISU']})
        >>> print(len(Isu_index))
        16
        >>> prot_coor.correct_cys_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['DISU']})
        >>> print(len(Isu_index))
        0
        """

        # FIND ISU res
        isu_index_list = self.get_index_selection({"res_name": ["DISU"]})
        if not isu_index_list:
            # print("Nothing to Fix")
            return self

        # Replace CYS resname :

        for atom_num in isu_index_list:
            self.atom_dict[atom_num]["res_name"] = "CYS"
            if self.atom_dict[atom_num]["name"] == "1CB":
                self.atom_dict[atom_num]["name"] = "CB"
            if self.atom_dict[atom_num]["name"] == "1SG":
                self.atom_dict[atom_num]["name"] = "SG"

        return self

    def water_to_ATOM(self):
        """ Change `HETATM` field of water to `ATOM`, as pdb2pqr only use ATOM field.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1dpx.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1dpx.pdb ,  1192 atoms found
        >>> hetatm_index = prot_coor.get_index_selection({'field':['HETATM']})
        >>> print(len(hetatm_index))
        179
        >>> prot_coor.water_to_ATOM() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> hetatm_index = prot_coor.get_index_selection({'field':['HETATM']})
        >>> print(len(hetatm_index))
        2
        >>> water_index = prot_coor.get_index_selection({'res_name':['HOH']})
        >>> print(len(water_index))
        177
        """

        # FIND Water res
        water_index_list = self.get_index_selection({'res_name':['HOH'], 'field':['HETATM']})
        if not water_index_list:
            return self
        else:
            self.change_index_pdb_field(water_index_list, {'field':'ATOM'})
            return self

    def correct_water_name(self):
        """ Correct the water resname from pdb2pqr

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1dpx.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1dpx.pdb ,  1192 atoms found
        >>> prot_coor.water_to_ATOM() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.write_pdb(os.path.join(TEST_OUT, '1dpx_water.pdb')) #doctest: +ELLIPSIS
        Succeed to save file .../1dpx_water.pdb
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_OUT, '1dpx_water.pdb'), os.path.join(TEST_OUT, '1dpx_water.pqr')) #doctest: +ELLIPSIS
        Succeed to read file .../1dpx_water.pdb ,  1192 atoms found
        Succeed to save file .../tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain .../tmp_pdb2pqr.pdb .../1dpx_water.pqr
        0
        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_OUT, '1dpx_water.pqr'), pqr_format = True) #doctest: +ELLIPSIS
        Succeed to read file .../1dpx_water.pqr ,  2491 atoms found
        >>> water_index = prot_coor.get_index_selection({'res_name':['TP3M'], 'name':['OH2']})
        >>> print(len(water_index))
        177
        """

        # FIND Water res
        water_index_list = self.get_index_selection({"res_name": ["TP3M"]})
        if not water_index_list:
            #print("Nothing no water to fix")
            return self

        # Replace SOL resname :

        for atom_num in water_index_list:
            self.atom_dict[atom_num]["res_name"] = "SOL"
            if self.atom_dict[atom_num]["name"] == "OH2":
                self.atom_dict[atom_num]["name"] = "OW"
            if self.atom_dict[atom_num]["name"] == "H1":
                self.atom_dict[atom_num]["name"] = "HW1"
            if self.atom_dict[atom_num]["name"] == "H2":
                self.atom_dict[atom_num]["name"] = "HW2"

        return self

    def insert_mol(self, pdb_out, out_folder, mol_chain, check_file_out=True):
        """
        Insert molecules defined by chain ID ``mol_chain`` in a water solvant.
        Check which water molecules are within ``cutoff_prot_off=12.0`` X
        and ``cutoff_prot_in=15.0`` X of protein and peptide C alpha atoms.
        Move the molecules to be inserted at the position of water molecules.
        Then delete all water molecules within ``cutoff_water_clash=1.2`` X of
        the inserted molecule atoms.

        :param pdb_out: name of output pdb file
        :type pdb_out: str

        :param out_folder: path of the ouput directory
        :type out_folder: str

        :param mol_chain: chain ID of the molecule to be inserted,
        :type mol_chain: str

        :param check_file_out: flag to check or not if file has already been created.
            If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        .. warning::
            self.atom_dict file must contain alredy a concatenated system with a
            ligand (chain: ``mol_chain``) and a solvated system.
        """

        # Create the out_folder:
        pdb_out = os.path.join(out_folder, pdb_out)
        os_command.create_dir(out_folder)

        # Parameters for molecule insertion:
        cutoff_water_clash = 1.2
        cutoff_prot_off = 12.0
        cutoff_prot_in = 15.0

        print("Insert mol in system")

        # Check if output files exist:
        if check_file_out and os.path.isfile(pdb_out):
            print("Insert Mol", pdb_out, "already exist")
            return None

        # Select protein, water and molecule atoms :
        prot_insert_CA = self.select_part_dict(selec_dict={'name': ['CA']})
        water = self.select_part_dict(selec_dict={'res_name': ['SOL']})
        water_O = self.select_part_dict(selec_dict={'res_name': ['SOL'], 'name': ['OW']})
        insert = self.select_part_dict(selec_dict={'chain': [mol_chain]})
        insert_ACE_C = self.select_part_dict(selec_dict={'chain': [mol_chain],
                                                         'name': ['C'],
                                                         'res_name': ['ACE']})

        mol_num = len(insert_ACE_C.atom_dict)
        res_insert_list = insert.get_attribute_selection(attribute='uniq_resid')
        # Need to sort the resid, to have consecutive residues
        res_insert_list.sort()
        #print("Residue list = ", res_insert_list)
        mol_len = int(len(res_insert_list) / mol_num)

        print("Insert {} mol of {:d} residues each".format(mol_num, mol_len))
        start_time = time.time()
        # Insert one molecule at a time:
        for i in range(mol_num):
            start_time = time.time()
            water_good_index = water_O.get_index_dist_between(prot_insert_CA,
                                                              cutoff_max=cutoff_prot_in,
                                                              cutoff_min=cutoff_prot_off)
            print('insert mol {:3}, water mol {:5}, time={:.2f}'.format(i+1, len(water_good_index), time.time() - start_time))
            insert_unique = insert.select_part_dict(selec_dict={'chain': [mol_chain],
                                                                'uniq_resid': res_insert_list[(mol_len * i):(mol_len * (i + 1))]})
            com_insert = insert_unique.center_of_mass()
            trans_vector = self.atom_dict[water_good_index[0]]['xyz'] - com_insert

            insert_unique.translate(trans_vector)

        # Delete water residues in which at leat one atom is close enough to peptides
        water_to_del_index = water.get_index_dist_between(insert, cutoff_max=cutoff_water_clash)
        water_res_to_del = water.get_attribute_selection(selec_dict={},
                                                         attribute='uniq_resid',
                                                         index_list=water_to_del_index)
        water_to_del_index = water.get_index_selection(selec_dict={'uniq_resid': water_res_to_del})
        print("Delete {} overlapping water atoms".format(len(water_to_del_index)))
        self.del_atom_index(index_list=water_to_del_index)

        self.write_pdb(pdb_out)
        return self

    def translate(self, vector):
        """ Translate all atoms of a coor object by a given ``vector``

        :param vector: 3d translation vector
        :type vector: list

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.01 y:0.45 z:8.57
        >>> prot_coor.translate(-com_1y0m)
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:-0.00 y:0.00 z:0.00

        """
        for atom_num, atom in self.atom_dict.items():
            atom['xyz'] += vector

    def center_of_mass(self, selec_dict={}):
        """ Compute the center of mass of a selection
        Avoid using atoms with 2 letters atom name like NA Cl ...

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={}

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.01 y:0.45 z:8.57

        .. warning::
            Atom name must start with its type letter (H, C, N, O, P, S).
        """

        com_array = np.zeros(3)
        mass_tot = 0

        atom_mass_dict = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16}

        for atom_num, atom in self.atom_dict.items():
            selected = True
            # print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("selection=",selection)
                # print("atom:",atom)
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                # print(atom)
                if atom['name'][0] in atom_mass_dict:
                    mass = atom_mass_dict[atom['name'][0]]
                    com_array += atom['xyz'] * mass
                    mass_tot += mass

        return com_array / mass_tot

    def get_index_dist_between(self, atom_sel_2, cutoff_min=0, cutoff_max=10):
        """ Check is distance between atoms of self.atom_dict is under cutoff
        with the atoms of group 1.
        Then return list of index of atoms of self.coor under cutoff ditance.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param cutoff_min: maximum distance cutoff
        :type cutoff_min: float, default=0.0

        :param cutoff_max: minimum distance cutoff
        :type cutoff_max: float, default=10.0
        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb')) #doctest: +ELLIPSIS
        Succeed to read file ...test/input/1y0m.pdb ,  648 atoms found
        >>> res_810 = prot_coor.select_part_dict({'res_num' : [810]})
        >>> close_r810 = prot_coor.get_index_dist_between(res_810, cutoff_min = 3, cutoff_max = 5)
        >>> print(len(close_r810))
        65

        """

        coor_array = np.array([atom['xyz'] for key, atom in self.atom_dict.items()])
        index_array = np.array([key for key, atom in self.atom_dict.items()])

        coor_array_2 = np.array([atom['xyz'] for key, atom in atom_sel_2.atom_dict.items()])
        
        #Compute distance matrix
        dist_mat = distance_matrix(coor_array, coor_array_2)

        # Compute index of matrix column under cutoff_max and over cutoff_min: 
        dist_mat_good =  np.where( (dist_mat.min(1) < cutoff_max) &  (dist_mat.min(1) > cutoff_min) )[0]
        
        return index_array[dist_mat_good]

    def dist_under_index(self, atom_sel_2, cutoff=10.0):
        """ Check is distance between atoms of self.coor is under cutoff with
        atoms of group 1.
        Then return list of index of atoms of self.coor under ctuoff ditance.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param cutoff: distance cutoff
        :type cutoff: float, default=10.0
        """

        coor_array = np.array([atom['xyz'] for key, atom in self.atom_dict.items()])
        index_array = np.array([key for key, atom in self.atom_dict.items()])

        coor_array_2 = np.array([atom['xyz'] for key, atom in atom_sel_2.atom_dict.items()])

        #Compute distance matrix
        dist_mat = distance_matrix(coor_array, coor_array_2)

        # Compute index of matrix column under cutoff_max and over cutoff_min: 
        dist_mat_good =  np.where( (dist_mat.min(1) < cutoff))[0]

        return index_array[dist_mat_good]

    def make_peptide(self, sequence, pdb_out, check_file_out=True):
        """
        Create a linear peptide structure.

        :param sequence: peptide sequence
        :type sequence: str

        :param pdb_out: name of output pdb file
        :type pdb_out: str

        :param check_file_out: flag to check or not if file has already been created.
            If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        """

        # Create and go in out_folder:
        # This is necessary for the topologie creation
        out_folder = os.path.dirname(pdb_out)
        # print(out_folder)
        os_command.create_dir(out_folder)


        print("-Make peptide: {}".format(sequence))

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            print("make_peptide", pdb_out, "already exist")
            return os.path.abspath(pdb_out)

        pep = Coor()
        seq = 'X' + sequence
        atom_num = 0
        uniq_resid = 0
        connect_dict = {}
        prev_res_name_index = {}

        # Initialize atom_dict:
        for res_name in seq:
            print("residue name:{}".format(res_name))
            res_name_index = {}

            for atom_name in AA_ATOM_DICT[res_name]:
                # print("\tatom name:{}".format(atom_name))

                if atom_num == 0:
                    xyz = np.zeros(3)
                else:
                    # Look for the previous bonded atom:
                    for dist in AA_BOND_DICT[res_name]:
                        if atom_name == dist[0]:
                            if dist[1][0] != '-':
                                connect_name = dist[1]
                                connect_index = res_name_index[connect_name]
                            else:
                                connect_name = dist[1][1:]
                                connect_index = prev_res_name_index[connect_name]
                            break
                        elif atom_name == dist[1]:
                            if dist[0][0] != '-':
                                connect_name = dist[0]
                                connect_index = res_name_index[connect_name]
                            else:
                                connect_name = dist[0][1:]
                                connect_index = prev_res_name_index[connect_name]
                            break
                    # print("{} connect to {} for {}".format(
                    #    atom_name, connect_name, res_name))
                    bond_len = Coor.find_dist(res_name, atom_name, connect_name)
                    # print("Bond : {}-{} = {} X".format(
                    #    atom_name, connect_name, bond_len))

                    if atom_num == 1:
                        xyz = pep.atom_dict[connect_index]['xyz'] + [bond_len, 0, 0]
                    connect_dict[atom_num] = connect_index

                    if atom_num > 1:
                        connect_2_index = connect_dict[connect_index]
                        connect_2_name = pep.atom_dict[connect_2_index]['name']
                        angle = Coor.find_angle(res_name, atom_name, connect_name, connect_2_name)
                        angle_rad = np.deg2rad(angle)
                        # print("Angle: {}-{}-{} = {}°".format(atom_name, connect_name, connect_2_name, angle))
                        if atom_num == 2:
                            xyz = pep.atom_dict[connect_index]['xyz'] + [
                                -bond_len * np.cos(np.deg2rad(angle)),
                                -bond_len * np.sin(np.deg2rad(angle)),
                                0]

                    if atom_num > 2:
                        if connect_2_index not in connect_dict:
                            xyz = pep.atom_dict[connect_index]['xyz'] + [-bond_len * np.cos(np.deg2rad(angle)),
                                                                         -bond_len * np.sin(np.deg2rad(angle)),
                                                                         0]
                        else:
                            connect_3_index = connect_dict[connect_2_index]
                            connect_3_name = pep.atom_dict[connect_3_index]['name']
                            dihe = Coor.find_dihe_angle(res_name, atom_name, connect_name, connect_2_name, connect_3_name)
                            dihe_rad = np.deg2rad(dihe)
                            # print("Dihedral Angle: {}-{}-{}-{} = {}°".format(atom_name, connect_name, connect_2_name, connect_3_name, dihe))
                            # From https://github.com/ben-albrecht/qcl/blob/master/qcl/ccdata_xyz.py#L208
                            vec_1 = pep.atom_dict[connect_index]['xyz'] - pep.atom_dict[connect_2_index]['xyz']
                            vec_2 = pep.atom_dict[connect_index]['xyz'] - pep.atom_dict[connect_3_index]['xyz']

                            vec_n = np.cross(vec_1, vec_2)
                            vec_nn = np.cross(vec_1, vec_n)

                            vec_n /= norm(vec_n)
                            vec_nn /= norm(vec_nn)

                            vec_n *= -sin(dihe_rad)
                            vec_nn *= cos(dihe_rad)

                            vec_3 = vec_n + vec_nn
                            vec_3 /= norm(vec_3)
                            vec_3 *= bond_len * sin(angle_rad)

                            vec_1 /= norm(vec_1)
                            vec_1 *= bond_len * cos(angle_rad)

                            xyz = pep.atom_dict[connect_index]['xyz'] + vec_3 - vec_1
                            # print(xyz)

                atom = {"field": 'ATOM',
                        "num": atom_num,
                        "name": atom_name,
                        "alter_loc": "",
                        "res_name": AA_1_TO_3_DICT[res_name],
                        "chain": "P",
                        "res_num": uniq_resid,
                        "uniq_resid": uniq_resid,
                        "insert_res": "",
                        "xyz": xyz,
                        "occ": 0.0,
                        "beta": 0.0}
                res_name_index[atom_name] = atom_num
                # print(atom)
                pep.atom_dict[atom_num] = atom
                atom_num += 1
            prev_res_name_index = res_name_index
            uniq_resid += 1

        pep.write_pdb(pdb_out)
        pdb_out = os.path.abspath(pdb_out)

        return pdb_out

    @staticmethod
    def find_dist(aa_name, name_a, name_b):
        for dist in DIST_DICT[aa_name]:
            if dist[:2] == [name_a, name_b] or dist[:2] == [name_b, name_a]:
                return dist[2]
        raise ValueError('Distance param {}-{} for {} not found !!'.format(name_a, name_b, aa_name))

    @staticmethod
    def find_angle(aa_name, name_a, name_b, name_c):
        for angle in ANGLE_DICT[aa_name]:
            if angle[:3] == [name_a, name_b, name_c] or angle[:3] == [name_c, name_b, name_a]:
                return angle[3]
        raise ValueError('Angle param {}-{}-{} for {} not found !!'.format(name_a, name_b, name_c, aa_name))

    @staticmethod
    def find_dihe_angle(aa_name, name_a, name_b, name_c, name_d):
        for angle in DIHE_DICT[aa_name]:
            if angle[:4] == [name_a, name_b, name_c, name_d] or angle[:4] == [name_d, name_c, name_b, name_a]:
                return angle[4]
        raise ValueError('Angle param {}-{}-{}-{} for {} not found !!'.format(name_a, name_b, name_c, name_d, aa_name))

    @staticmethod
    def atom_dist(atom_1, atom_2):
        """Compute the distance between 2 atoms.

        :param atom_1: atom dictionnary
        :type atom_1: dict

        :param atom_2: atom dictionnary
        :type atom_2: dict

        :return: distance
        :rtype: float

        :Example:

        >>> atom_1 = {'xyz': np.array([0.0, 0.0, 0.0])}
        >>> atom_2 = {'xyz': np.array([0.0, 1.0, 0.0])}
        >>> atom_3 = {'xyz': np.array([1.0, 1.0, 1.0])}
        >>> Coor.atom_dist(atom_1, atom_2)
        1.0
        >>> Coor.atom_dist(atom_1, atom_3)
        1.7320508075688772
        """

        distance = np.linalg.norm(atom_1['xyz'] - atom_2['xyz'])
        return distance

    @staticmethod
    def concat_pdb(*pdb_in_files, pdb_out):
        """Concat a list of pdb files in one.

        :param pdb_in_files: list of pdb files
        :type pdb_in_files: list

        :param pdb_out: atom dictionnary
        :type pdb_out: dict

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> Coor.concat_pdb(os.path.join(TEST_PATH, '1y0m.pdb'),
        ...                 os.path.join(TEST_PATH, '1rxz.pdb'),
        ...                 pdb_out = os.path.join(TEST_OUT, 'tmp_2.pdb')) #doctest: +ELLIPSIS
        Succeed to save concat file: .../tmp_2.pdb
        """

        if os_command.check_file_and_create_path(pdb_out):
            print("File " + pdb_out + " already exist")
            return

        filout = open(pdb_out, 'w')
        count = 0

        for pdb_in in pdb_in_files:
            #print("Concat :", os.path.relpath(pdb_in))
            with open(pdb_in) as pdbfile:
                for line in pdbfile:
                    if (count == 0 and line[:6] == "CRYST1") or line[:4] == 'ATOM' or line[:6] == "HETATM":
                        filout.write(line)
            count += 1

        filout.close()
        print("Succeed to save concat file: ", pdb_out)


if __name__ == "__main__":

    import doctest
    import shutil

    TEST_DIR = 'gromacs_py_test_out'
    TEST_OUT = os.path.join(TEST_DIR, 'pdb_manip_test')

    def getfixture(*args):
        return TEST_OUT

    #sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    #print(os.path.abspath(os.path.join(os.path.dirname(__file__))))

    print("-Test pdb_manip module:")
    print("pdb_manip:\t", doctest.testmod())

    # Erase all test files
    shutil.rmtree(TEST_DIR, ignore_errors=True)
