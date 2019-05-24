Installation
=======================================

Conda installation
---------------------------------------

Quick Start
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then add several conda channels for dependencies:

.. code-block:: bash

	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda config --add channels samuel.murail

Finally install `gromacs_py` in your `gromacs_py_env` environment:

.. code-block:: bash

	conda install gromacs_py

In a new environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First create a conda environment named `gromacs_py_env` (choose yout env name). `gromacs_py` support python version 3.5 to 3.7.

.. code-block:: bash

	conda create --yes -n gromacs_py_env python=3.7

Then add several conda channels for dependencies:
	- gromacs 2019 from `bioconda`
	- gromacs_py from my channel (`samuel.murail`)
	- htmd-pdb2pqr from `conda-forge`


.. code-block:: bash

	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda config --add channels samuel.murail

Finally install `gromacs_py` in your `gromacs_py_env` environment:

.. code-block:: bash

	conda install --yes -n gromacs_py_env gromacs_py

You need to activate the environment to be able to use `gromacs_y`, it has to be done in every shell in which you need `gromacs_py` :

.. code-block:: bash

	source activate gromacs_py_env

Finally test the installation using `pytest`:

.. code-block:: bash

	(gromacs_py_env) $ pip install pytest
	(gromacs_py_env) $ pytest --pyargs gromacs_py.gromacs --doctest-modules


Without Conda
---------------------------------------



Get the gromacs_py library from `github`_.

.. code-block:: bash

	git clone https://github.com/samuelmurail/gromacs_py.git
	./setup.py install --user

.. _github: https://github.com/samuelmurail/gromacs_py

Prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. python 3 libraries:  
	* numpy
	* scipy
	* pandas
	* matplotlib
	* Sphinx and sphinx-argparse (only for building documentation)

2. `pdb2pqr`_:

.. code-block:: bash

	git clone https://github.com/Electrostatics/apbs-pdb2pqr.git --branch master --depth=1
	cd apbs-pdb2pqr/pdb2pqr/
	python scons/scons.py install --prefix=$HOME

3. `Gromacs`_

Get source code from `gromacs website`__ and follow the following command for a quick and dirty install (for more details see `gromacs 2019 install guide`__)

In my case I add to change few options to ``cmake``:

	* ``-DCMAKE_C_COMPILER=gcc-6``, as gcc versions later than 6 are not supported.
	* ``-DGMX_GPU=on`` to use GPU acceleration
	* ``-DCMAKE_INSTALL_PREFIX=../../local-gromacs-2019.2/`` to install gromacs in a non-standard location

.. code-block:: bash

	tar -xfz gromacs-2019.2.tar.gz
	cd gromacs-2019.2
	mkdir build
	cd build
	cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=gcc-6 -DGMX_GPU=on -DCMAKE_INSTALL_PREFIX=../../local-gromacs-2019.2/ 

	# the option -j 4 allow using 4 processor for compilation
	make -j 4
	make check -j 4
	make install -j 4
	
	source ../../local-gromacs-2019.2/bin/GMXRC


.. _pdb2pqr: http://www.poissonboltzmann.org/
.. _Gromacs: http://www.gromacs.org/
__ http://manual.gromacs.org/documentation/
__ http://manual.gromacs.org/documentation/2019/install-guide/index.html

Installing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Need to add path of gmx and pdb2pqr to the environment variable ``$PATH``.
Add in your ~/.bashrc :

.. code-block:: bash

	# Add gromacs 'gmx' path:
	export PATH='*path_to_gromacs*/bin/':$PATH
	# Add pdb2pqr 'pdb2pqr.py' path:
	export PATH='*path_to_apbs-pdb2pqr/pdb2pqr/':$PATH


Make the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Need `sphinx`_ installed with the argparse sphinx module:

.. code-block:: bash

	pip3 install Sphinx --user
	pip3 install sphinx-argparse --user

You can then build the documentation either in html format or pdf.

.. code-block:: bash

	cd gromacs_py/doc
	# For html documentation:
	sphinx-build -b html . _build
	# For pdf documentation:
	sphinx-build -M latexpdf . _build/

.. _sphinx: http://www.sphinx-doc.org

Test installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Launch test with `doctest`_, will check that module’s docstrings are up-to-date by verifying that all interactive examples still work as documented.

.. code-block:: bash

	$ ./test_gromacs_py.py
	tools.os_command:  	 TestResults(failed=0, attempted=19)
	tools.pdb_manip:	 TestResults(failed=0, attempted=127)
	tools.pdb2pqr:  	 TestResults(failed=0, attempted=11)
	gromacs.gmx5:    	 TestResults(failed=0, attempted=52)

.. _doctest: https://docs.python.org/3/library/doctest.html