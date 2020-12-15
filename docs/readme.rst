.. gromacs_py documentation master file, created by
   sphinx-quickstart on Fri Jun  1 16:36:00 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. include:: ../README.rst

Main features:
---------------------------------------

* Python Scriptable simulation:
   * Topologie creation
   * Solvation
   * Ion insertion
   * Energy minimisation
   * Equilibration with different position restraints
   * Production

* Topologie manipulation starting from a raw ``PDB``:
   * Amino acid protonation and pKa calculation using `apbs/pdb2pqr <http://www.poissonboltzmann.org/>`_
   * Position constraints file ``.itp`` creation
   * Cyclic peptide topologie
   * Cystein bond topologie modification
   * ligand topologie using `ambertools` and `acpype`

* Advanced simulation tools:
   * Monitor a simulation while running
   * Free Energy calculations
   * Interrupt a simulation if a criterion is met (Not implemented yet)


Compatibility
---------------------------------------

* Supported Gromacs versions:
   * 2020
   * 2019*
   * 2018*
   * 2017
   * 2016
   * 5.1
   * 5.0

* Supported Python versions:
   * 3.8*
   * 3.7*
   * 3.6*

* Supported OS:
   * osx*
   * linux*

**\*** tested after each code submission.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
