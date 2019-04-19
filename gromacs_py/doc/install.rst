
Getting Started
=======================================

Get the gromacs_py from `github`_.

.. code-block:: bash

	git clone https://github.com/samuelmurail/gromacs_py.git
	./setup.py install

.. _github: https://github.com/samuelmurail/gromacs_py

Prerequisites
=======================================

1. python 3 librairies:  
	* sys, os, shutil, glob, argparse, subprocess, operator

2. `pdb2pqr`_:

.. code-block:: bash

	git clone https://github.com/Electrostatics/apbs-pdb2pqr.git --branch master --depth=1
	cd apbs-pdb2pqr/pdb2pqr/
	python scons/scons.py install

3.  `Gromacs`_

.. _pdb2pqr: http://www.poissonboltzmann.org/
.. _Gromacs: http://www.gromacs.org/


Installing
=======================================

Need to add path of gmx and pdb2pqr to the environment variable ``$PATH``.
Add in your ~/.bashrc :

.. code-block:: bash

	# Add gromacs 'gmx' path:
	export PATH='*path_to_gromacs*/bin/':$PATH
	# Add pdb2pqr 'pdb2pqr.py' path:
	export PATH='*path_to_apbs-pdb2pqr/pdb2pqr/':$PATH


Make the documentation
=======================================

Need sphinx install with argparse sphinx module.

.. code-block:: bash

	cd gromacs_py/doc
	# For html documentation:
	sphinx-build -b html . _build
	# For pdf documentation:
	sphinx-build -M latexpdf . _build/

Test installation
=======================================

Launch test with `doctest`__, will check that module’s docstrings are up-to-date by verifying that all interactive examples still work as documented.

.. code-block:: bash

	$ ./test_gromacs_py.py
	tools.os_command:  	 TestResults(failed=0, attempted=19)
	tools.pdb_manip:	 TestResults(failed=0, attempted=127)
	tools.pdb2pqr:  	 TestResults(failed=0, attempted=11)
	gromacs.gmx5:    	 TestResults(failed=0, attempted=52)

__  https://docs.python.org/3/library/doctest.html