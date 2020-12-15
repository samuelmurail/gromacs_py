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


**Gromacs_py** is a Python library allowing a simplified use of the Gromacs MD simulation software. **Gromacs_py** can build a system topologie based on a pdb file, create the simulation system (pbc box, adding water and ions) and run minimisation, equilibration and production runs.
One of the main objective of the **Gromacs_py** wrapper is to automatize routine operations for MD simulation of multiple systems.

**Gromacs_py** is under active development using continuous integration with `Travis Cl <https://travis-ci.org/samuelmurail/gromacs_py>`_. 

* Online Documentation:
   https://gromacs-py.readthedocs.io

* Source code repository:
   https://github.com/samuelmurail/gromacs_py

Quick install
---------------------------------------

The latest release can be installed via `pip` or `conda`.

Pip
***************************************

If Gromacs (version >= 5.1) is already install, then you need to install the **Gromacs_py** library using `pypi <https://pypi.org/project/gromacs-py/>`_, and add the Gromacs ``gmx`` command in the environmnent variable ``$PATH``:

.. code-block:: bash

   pip install gromacs_py

   # Add gromacs 'gmx' path:
   export PATH='*path_to_gromacs*/bin/':$PATH

Conda
***************************************

If you don't need a GPU compiled version of Gromacs you can use directly the **Gromacs_py** `conda package <https://anaconda.org/bioconda/gromacs_py>`_ to install both Gromacs software and **Gromacs_py** library:

.. code-block:: bash

   conda install -c bioconda gromacs_py

Authors
---------------------------------------

* `Samuel Murail <https://samuelmurail.github.io/PersonalPage/>`_, Associate Professor - `Université Paris Diderot <https://www.univ-paris-diderot.fr>`_, `CMPLI <http://bfa.univ-paris-diderot.fr/equipe-8/>`_.

See also the list of `contributors <https://github.com/samuelmurail/gromacs_py/contributors>`_ who participated in this project.

License
---------------------------------------

This project is licensed under the GNU General Public License v2.0 - see the ``LICENSE`` file for details.
