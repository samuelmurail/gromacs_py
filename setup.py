#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Create the pipy:
# python setup.py sdist
# twine check dist/*
# twine upload dist/*

# Create conda env:
# conda update conda
# conda env create --name test -f=.conda.yml
# conda activate test
# conda deactivate

from setuptools import setup, find_packages


setup(
    name='gromacs_py',
    version='1.1.1',
    author='Samuel Murail',
    author_email='samuel.murail@u-paris.fr',
    packages=find_packages(),
    scripts=['bin/extend.py',
             'bin/create_peptide.py',
             'bin/create_top.py',
             'bin/prepare_prot.py',
             'bin/insert_mol_no_vmd.py',
             'bin/minimize_pdb.py',
             'bin/solvate_ions.py',
             'bin/prepare_prot_topo_edit.py',
             'bin/production.py',
             'bin/test_gromacs_py.py',
             'bin/top_em_equi_3_step_prod.py'],
    description='Gromacs_py is a python library allowing a simplified use of the gromacs MD simulation software.',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type='text/x-rst',
    url='https://github.com/samuelmurail/gromacs_py',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.5',
    install_requires=[
        "pandas",
        "matplotlib",
        "os_command_py",
        "pdb2pqr_htmd_propka30",
        "pdb_manip_py",
    ],
    package_data={
        '': ['docs/*rst',
             'docs/*py',
             'docs/notebook/00_basic_example.ipynb'],
        'gromacs_py': ['template/*mdp',
                       'template/charmm36-jul2017.ff/*',
                       'test_files/*pdb',
                       'test_files/*xvg'],
    },
    project_urls={
        'Bug Reports': 'https://github.com/samuelmurail/gromacs_py/issues',
        'Funding': 'https://www.impots.gouv.fr/portail/',
        'Source': 'https://github.com/samuelmurail/gromacs_py',
    },
)