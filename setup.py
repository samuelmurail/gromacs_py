#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import gromacs_py

setup(
    name='gromacs_py',
    version=gromacs_py.__version__,
    packages=find_packages(),
    description='Gromacs_py is a python library allowing a simplified use of the gromacs MD simulation software.',
    long_description=open('README.rst', encoding='utf-8').read(),
    author='Samuel Murail',
    author_email='samuel.murail at univ-paris-diderot.fr',
    url='https://github.com/samuelmurail/gromacs_py',
    classifiers=[
        "Development Status :: stable",
        "License :: GPL-2.0",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.5",
        "Topic :: MD simulation",
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
