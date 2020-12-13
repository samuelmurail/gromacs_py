.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1455734.svg
   :target: https://doi.org/10.5281/zenodo.1455734

.. image:: https://readthedocs.org/projects/gromacs-py/badge/?version=latest
   :target: https://gromacs-py.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.org/samuelmurail/gromacs_py.svg?branch=master
   :target: https://travis-ci.org/samuelmurail/gromacs_py

.. image:: https://codecov.io/gh/samuelmurail/gromacs_py/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/samuelmurail/gromacs_py

.. image:: https://anaconda.org/bioconda/gromacs_py/badges/version.svg
   :target: https://anaconda.org/bioconda/gromacs_py


.. image:: https://badge.fury.io/py/gromacs-py.svg
   :target: https://badge.fury.io/py/gromacs-py

Gromacs_py
=======================================


Gromacs_py is a python library allowing a simplified use of the gromacs MD simulation software. Gromacs_py can build topologie based on a pdb file, create the simulation system (box, add water and ions) and run minimisation, equilibration and production run.
One of the main objective of the gromacs_py wrapper is to automatize routine operations for MD simulation of multiple systems.

* Online Documentation:
   https://gromacs-py.readthedocs.io

* Source code repository:
   https://github.com/samuelmurail/gromacs_py

Main features:
---------------------------------------

* Python Scriptable simulation:
   - Topologie creation
   - Solvation
   - Ion insertion
   - Energy minimisation
   - Equilibration with different position restraints
   - Production

* Topologie manipulation starting from a raw ``PDB``:
   - Amino acid protonation and pKa calculation using `apbs/pdb2pqr <http://www.poissonboltzmann.org/>`_
   - Position constraints file ``.itp`` creation
   - Cyclic peptide topologie
   - Cystein bond topologie modification
   - ligand topologie using `ambertools` and `acpype`

* Advanced simulation tools:
   - Monitor a simulation while running
   - Free Energy calculations
   - Interrupt a simulation if a criterion is met (Not implemented yet)


Compatibility
---------------------------------------

* Supported Gromacs versions:
   - 2020
   - 2019*
   - 2018*
   - 2017
   - 2016
   - 5.1
   - 5.0

* Supported Python versions:
   - 3.8*
   - 3.7*
   - 3.6*
   - 3.5*

* Supported OS:
   - osx*
   - linux*

**\*** tested after each code submission.

Quick install
---------------------------------------

With pip
***************************************

If gromacs (version >= 5.1) is already install, then install you need to install the `gromacs_py` library, and add the gromacs `gmx` command in the environmnent variable `$PATH`:

.. code-block:: bash

   pip install gromacs_py

   # Add gromacs 'gmx' path:
   export PATH='*path_to_gromacs*/bin/':$PATH

With conda
***************************************

If you don't need a GPU version you can use directly the `gromacs_py` conda package:

.. code-block:: bash

   conda install gromacs_py

Authors
---------------------------------------

* `Samuel Murail <https://samuelmurail.github.io/PersonalPage/>`_, Associate Professor - `Université Paris Diderot <https://www.univ-paris-diderot.fr>`_, `CMPLI <http://bfa.univ-paris-diderot.fr/equipe-8/>`_.

See also the list of `contributors <https://github.com/samuelmurail/gromacs_py/contributors>`_ who participated in this project.

License
---------------------------------------

This project is licensed under the GNU General Public License v2.0 - see the ``LICENSE`` file for details.
