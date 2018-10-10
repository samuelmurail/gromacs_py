# Gromacs_lib

Gromacs_lib is a python library allowing a simplified use of the gromacs MD simulation software. Gromacs_lib can build topologie based on a pdb file, create the simulation system (box, add water and ions) and run minimisation, equilibration and production run.

## Getting Started
Get the gromacs_py from the [RPBS](https://gitlab.rpbs.univ-paris-diderot.fr) gitlab.
```bash
git clone https://github.com/samuelmurail/gromacs_py.git
```


### Prerequisites

1.  python librairies:  
   - sys, os, shutil, glob, argparse, subprocess

2.  [pdb2pqr](http://www.poissonboltzmann.org/):
```bash
git clone https://github.com/Electrostatics/apbs-pdb2pqr.git --branch master --depth=1
cd apbs-pdb2pqr/pdb2pqr/
python scons/scons.py install
```
3.  [Gromacs](http://www.gromacs.org/)

4.  [VMD](http://www.ks.uiuc.edu/Research/vmd/)  
   - Only used to insert peptide/ligands is the simulation box.

5.  [Pymol](https://pymol.org/2/)  
   - Pymol is only used for *de novo* peptide creation.


### Installing

Explain how to setup gromacs lib **TO DO**
```
Give the example
```

## Running the tests

Explain how to run the automated tests for this system
```bash
cd gromacs_py/test
test_gromacs5.py
test_pdb.py
test_cyclic_peptide.py
```

## Usage: **TO DO** 
#### basic:
1.  create_top.py
2.  solvate_ions.py
3.  minimize_pdb.py
4.  equi_3_step.py
5.  production.py
6.  extend.py
7.  mini_equi_3_step.py (*Remove ?*)  

#### optional:
1. concat_pdb.py
2. insert_mol.py
3. create_peptide.py
4. sim_pep_prot.py



## Authors

* **Samuel Murail** - Université Paris Diderot, [MTi lab](http://www.mti.univ-paris-diderot.fr/) 

**TO DO** See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

**TO DO** This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

**TO DO**