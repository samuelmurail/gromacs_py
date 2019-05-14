#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import gromacs_py

setup(
    name='gromacs_py',
    version=gromacs_py.__version__,
    packages=find_packages(),
    description='Gromacs wrapper library',
    long_description=open('README.rst', encoding='utf-8').read(),
    author='Samuel Murail',
    author_email='samuel.murail at univ-paris-diderot.fr',
    url='https://github.com/samuelmurail/gromacs_py',
    classifiers=[
        "Programming Language :: Python3",
        "Development Status :: alpha",
        "License :: To determine",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
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
        "pandoc"
    ],
    include_package_data=True,
)
