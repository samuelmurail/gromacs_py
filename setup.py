#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Create the pipy:
# python setup.py sdist
# twine check dist/*
# twine upload dist/*

from setuptools import setup, find_packages

import gromacs_py

setup(
    name='gromacs_py',
    version=gromacs_py.__version__,
    packages=find_packages(),
    description='Gromacs_py is a python library allowing a simplified use of the gromacs MD simulation software.',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type='text/x-rst',
    author='Samuel Murail',
    author_email='samuel.murail@univ-paris-diderot.fr',
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
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "sphinx_rtd_theme",
        "sphinx-argparse",
        "nbsphinx",
    ],
    include_package_data=True,
)
