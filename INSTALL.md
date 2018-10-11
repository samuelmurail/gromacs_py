
## Getting Started

Get the gromacs_py from [github](https://github.com/samuelmurail/gromacs_py).

```bash
git clone https://github.com/samuelmurail/gromacs_py.git
./setup.py install
```



## Prerequisites

1. python 3 librairies:  
	* sys, os, shutil, glob, argparse, subprocess, operator

2. [pdb2pqr](http://www.poissonboltzmann.org/):
```bash
git clone https://github.com/Electrostatics/apbs-pdb2pqr.git --branch master --depth=1
cd apbs-pdb2pqr/pdb2pqr/
python scons/scons.py install
```

3.  [Gromacs](http://www.gromacs.org/)

4.  [VMD](http://www.ks.uiuc.edu/Research/vmd/)
	* Only used to insert peptide/ligands is the simulation box.

5.  [Pymol](https://pymol.org/2/)
	* Pymol is only used for *de novo* peptide creation.



## Installing

Need to add path of gmx, pdb2pqr.py, vmd and pymol to the environment variable ``$PATH``.
Add in your ~/.bashrc :

```bash
# Add gromacs 'gmx' path:
export PATH='*path_to_gromacs*/bin/':$PATH
# Add VMD 'vmd_MACOSXX86' or 'vmd' path:
export PATH='/Applications/VMD\ 1.9.3.app/Contents/vmd/':$PATH
# Add PyMol 'pymol' path:
export PATH='/Applications/PyMOL.app/Contents/bin/':$PATH
# Add pdb2pqr 'pdb2pqr.py' path:
export PATH='*path_to_apbs-pdb2pqr/pdb2pqr/':$PATH
```


## Make the documentation

Need sphinx install with m2r module.

```bash
cd gromacs_py/doc
# For html documentation:
sphinx-build -b html . _build
# For pdf documentation:
sphinx-build -M latexpdf . _build/
```

## Test installation

Launch test with [doctest](https://docs.python.org/3/library/doctest.html), will check that module’s docstrings are up-to-date by verifying that all interactive examples still work as documented.

```bash
$ ./test_gromacs_py.py
tools.pdb_manip:	 TestResults(failed=0, attempted=88)
tools.pdb2pqr:  	 TestResults(failed=0, attempted=11)
gromacs.gmx5:    	 TestResults(failed=0, attempted=52)
```


