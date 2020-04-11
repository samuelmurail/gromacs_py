

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1455734.svg
   :target: https://doi.org/10.5281/zenodo.1455734


.. image:: https://readthedocs.org/projects/gromacs-py/badge/?version=latest
   :target: https://gromacs-py.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. image:: https://travis-ci.org/samuelmurail/gromacs_py.svg?branch=master
    :target: https://travis-ci.org/samuelmurail/gromacs_py

.. image:: https://codecov.io/gh/samuelmurail/gromacs_py/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/samuelmurail/gromacs_py

.. image:: https://anaconda.org/samuel.murail/gromacs_py/badges/version.svg
  :target: https://anaconda.org/samuel.murail/gromacs_py

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
   - Equilibration with position restraints
   - Production
   - **GPU** acceleration

* Topologie manipulation starting from a raw ``PDB``:
   - Amino acid protonation and pKa calculation using `apbs/pdb2pqr <http://www.poissonboltzmann.org/>`_
   - Position constraints file ``.itp`` creation
   - Cyclic petide topologie

* Coordinate manipulation:
   - Changing atom names, chain, coordinates, ...
   - Insertion of *N* copy of a molecule in a system for *flooding* simulation
   - Linear peptide creation

* Advanced simulation tools:
   - Monitor a simulation while running
   - Interrupt a simulation if a criterion is met


Compatibility
---------------------------------------

* Supported Gromacs versions:
   - 2019*
   - 2018*
   - 2017
   - 2016
   - 5.1
   - 5.0

* Supported Python versions:
   - 3.7*
   - 3.6*
   - 3.5*

* Supported OS:
   - osx*
   - linux*

**\*** tested after each code submission.

Quick install
---------------------------------------

Add several conda channels for dependencies:

.. code-block:: bash

   conda config --add channels conda-forge
   conda config --add channels bioconda
   conda config --add channels samuel.murail

Then install `gromacs_py`:

.. code-block:: bash

   conda install gromacs_py

see `Installation <https://gromacs-py.readthedocs.io/en/latest/install.html>`_ for details.

.. _github: https://github.com/samuelmurail/gromacs_py

Tutorial
---------------------------------------

Here is an example of a short simulation of the SH3 domain of phospholipase Cγ1.
Seven successive steps are used:

1. Topologie creation using ``create_top.py``.
2. Minimisation of the structure using ``minimize_pdb.py``.
3. Solvation of the system using ``solvate_ions.py``.
4. Minimisation of the system using ``minimize_pdb.py``.
5. Equilibration of the system using ``equi_3_step.py``.
6. Production run using ``production.py``.
7. Extension of the production run using ``extend.py``.

.. code-block:: bash

   # Create topologie
   gromacs_py/create_top.py -f gromacs_py/test/input/1y0m.pdb  -o tmp/1y0m/top -vsite

   # Minimize the protein structure
   gromacs_py/minimize_pdb.py -f tmp/1y0m/top/1y0m_pdb2gmx_box.pdb -p tmp/1y0m/top/1y0m_pdb2gmx.top -o tmp/1y0m/em/  -n em_1y0m -nt 2

   # Add water and ions
   gromacs_py/solvate_ions.py -f tmp/1y0m/em/em_1y0m_compact.pdb -p tmp/1y0m/top/1y0m_pdb2gmx.top -o tmp/1y0m_water_ions/top/  -n 1y0m_water_ions

   # Minimize the system
   gromacs_py/minimize_pdb.py -f tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.gro -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/em/  -n em_1y0m

   # Do three small equilibrations with postion contraints on heavy atoms (first), Carbon alpha (second) and low constraint on Carbon alpha (third)
   gromacs_py/equi_3_step.py -f tmp/1y0m_water_ions/em/em_1y0m_compact.pdb -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/  -n 1y0m -HA_time 0.1 -CA_time 0.1 -CA_LOW_time 0.1

   # Small production run of 0.1 ns
   gromacs_py/production.py -f tmp/1y0m_water_ions/02_equi_CA_LOW/equi_CA_LOW_1y0m.gro -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/03_prod -n 1y0m -time 0.1

   # Extension of the simulation
   gromacs_py/extend.py -s tmp/1y0m_water_ions/03_prod/prod_1y0m.tpr -time 0.2

   # Remove simulation files
   rm -r ./tmp

Or simply use one command to do all previous commands:

.. code-block:: bash

   gromacs_py/top_em_equi_3_step_prod.py -f gromacs_py/test/input/1y0m.pdb -o tmp/1y0m -vsite -HA_time 0.1 -CA_time 0.1 -CA_LOW_time 0.1 -prod_time 0.3

Authors
---------------------------------------

* `Samuel Murail <https://samuelmurail.github.io/PersonalPage/>`_, Associate Professor - `Université Paris Diderot <https://www.univ-paris-diderot.fr>`_, `CMPLI <http://bfa.univ-paris-diderot.fr/equipe-8/>`_.

See also the list of `contributors <https://github.com/samuelmurail/gromacs_py/contributors>`_ who participated in this project.

License
---------------------------------------

This project is licensed under the GNU General Public License v2.0 - see the ``LICENSE`` file for details.
