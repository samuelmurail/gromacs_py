.. highlight:: shell

============
Installation
============

1. Get sources from the `GithubRepo`_
--------------------------------------

The sources for Gromacs_py can be downloaded from the `GithubRepo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/samuelmurail/gromacs_py

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/samuelmurail/gromacs_py/tarball/master

Once you have a copy of the source, switch to the ``gromacs_py`` directory.

.. code-block:: console

    $ cd gromacs_py


.. _GithubRepo: https://github.com/samuelmurail/gromacs_py
.. _tarball: https://github.com/samuelmurail/gromacs_py/tarball/master


2. Create Conda Environment
---------------------------

You need to create a conda environment to be able to use:  

* `Gromacs`_
* `Rdkit`_ Used for ligand parametrization, convert SMILE to pdb.
* `Antechamber`_ Amber tools for ligand parametrization.
* `Acpype`_ a python tool to use antechamber.
* `Apbs Pdb2pqr`_ Protein protonation calculation.

.. _Gromacs: http://www.gromacs.org/
.. _Rdkit: https://www.rdkit.org/
.. _Antechamber: http://ambermd.org/antechamber/
.. _Acpype: https://github.com/alanwilter/acpype
.. _Apbs Pdb2pqr: https://www.poissonboltzmann.org/


Use `conda en create` to create it using the ``.conda.yml`` file. You can overide the environmnent name using the option ``--name YOUR_NAME``.

.. code-block:: console

    $ conda env create -f .conda.yml

If you plan to use ``gromacs_py`` in jupyter notebook, you should try the ``jupyter`` version:

.. code-block:: console

    $ conda env create -f .conda_jupyter.yml


This will create an environmnet called ``gromacs_py`` (or the name you defined). You will then, need to activate the environmnent:

.. code-block:: console

    $ conda activate gromacs_py

.. note::
	If you want to install yourself Gromacs to be able to use the GPU acceleration, you can use the ``.conda_no_gromacs.yml`` or ``.conda_jupyter_no_gromacs.yml``:

	.. code-block:: console

	    $ conda env create -f .conda_no_gromacs.yml
	    $ conda activate gromacs_py

3. Install gromacs_py
---------------------

Once you have a copy of the source and have create a conda encironment,
you can install it with:

.. code-block:: console

    $ python setup.py install



4. Test Installation
--------------------

To test the installation, simply use ``pytest``:

.. code-block:: bash

	$ pytest
	=========================== test session starts ========================
	platform linux -- Python 3.8.2, pytest-5.4.2, py-1.9.0, pluggy-0.13.1
	rootdir: /home/murail/Documents/Code/gromacs_py, inifile: pytest.ini
	plugins: cov-2.10.1
	collected 30 items

	gromacs_py/gmx.py .............                                   [ 43%]
	gromacs_py/test/test_FreeEner.py ......                           [ 63%]
	gromacs_py/test/test_GmxSys.py ..                                 [ 70%]
	gromacs_py/tools/ambertools.py ....                               [ 83%]
	gromacs_py/tools/monitor.py .....                                 [100%]

	======================= 30 passed in 236.83s (0:03:56) =================


Conda installation
---------------------------------------

If you don't need a GPU compiled version of Gromacs you can use directly the **Gromacs_py** `conda package <https://anaconda.org/bioconda/gromacs_py>`_ to install both Gromacs software and **Gromacs_py** library:

.. code-block:: bash

   conda install -c bioconda gromacs_py


Pypi (Deprecated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If gromacs (version >= 5.1) is already install, then install you need to install the `gromacs_py` library, and add the gromacs `gmx` command in the environmnent variable `$PATH`:

.. code-block:: bash

	pip install gromacs_py

	# Add gromacs 'gmx' path:
	export PATH='*path_to_gromacs*/bin/':$PATH


Without Conda 
---------------------------------------

Get the gromacs_py library from `github`_.

.. code-block:: bash

	git clone https://github.com/samuelmurail/gromacs_py.git
	./setup.py install --user

	# Add gromacs 'gmx' path:
	export PATH='*path_to_gromacs*/bin/':$PATH

.. _github: https://github.com/samuelmurail/gromacs_py

Prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. python 3 libraries installed when you launch the pip command:  
	* numpy
	* scipy
	* pandas
	* matplotlib
	* Sphinx and sphinx-argparse (only for building documentation)
	* `Os_Command_py`_
	* `PDB_Manip_py`_
	* `PDB2PQR`_ using the package `pdb2pqr_htmd_propka30`_ a python 3 version developped by `tonigi`_ and adapted to use successfully propka3.0.

.. _Os_Command_py: https://github.com/samuelmurail/os_command_py
.. _PDB_Manip_py: https://github.com/samuelmurail/pdb_manip_py
.. _PDB2PQR: http://www.poissonboltzmann.org/
.. _pdb2pqr_htmd_propka30: https://github.com/samuelmurail/apbs-pdb2pqr/tree/htmd-fixups
.. _tonigi: https://github.com/tonigi/apbs-pdb2pqr

2. `Gromacs`_

Get source code from `gromacs website`__ and follow the following command for a quick and dirty install (for more details see `gromacs 2019 install guide`__)

In my case I add to change few options to ``cmake``:

	* ``-DCMAKE_C_COMPILER=gcc-8``, as gcc versions later than 6 are not supported.
	* ``-DGMX_GPU=on`` to use GPU acceleration
	* ``-DCMAKE_INSTALL_PREFIX=../../local-gromacs-2019.2/`` to install gromacs in a non-standard location

.. code-block:: bash

	# Specify the version:
	version="2021.5"
	# To modify:
	dir_install="/home/murail/Documents/Software/local-gromacs-${version}"

	wget https://ftp.gromacs.org/gromacs/gromacs-${version}.tar.gz
	tar -xvf gromacs-${version}.tar.gz
	cd gromacs-${version}
	mkdir build
	cd build
	# In my case I needed to define ggc-8 because gromacs doesn't accept superior versions
	cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=${dir_install} -DCMAKE_C_COMPILER=gcc-8
	make
	make check
	make install
	source ${dir_install}/bin/GMXRC
	echo "export PATH=${dir_install}/bin/:\$PATH" >> ~/.bashrc


.. _Gromacs: http://www.gromacs.org/
__ http://manual.gromacs.org/documentation/
__ http://manual.gromacs.org/documentation/2019/install-guide/index.html

3. `Ambertools`_

The easiest way is to use the conda package. However is can also be installed from source.

.. _Ambertools: https://ambermd.org/AmberTools.php
