# Gromacs_py

Gromacs_py is a python library allowing a simplified use of the gromacs MD simulation software. Gromacs_py can build topologie based on a pdb file, create the simulation system (box, add water and ions) and run minimisation, equilibration and production run.

## Quick install

Get the gromacs_py from [github](https://gitlab.rpbs.univ-paris-diderot.fr).

```bash
	git clone git@gitlab.rpbs.univ-paris-diderot.fr:murail/gromacs_py.git
	./setup.py install
```


see the [INSTALL.md](INSTALL.md) file for details

## Usage:
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

* **Samuel Murail** - [Université Paris Diderot](https://www.univ-paris-diderot.fr), [MTi lab](http://www.mti.univ-paris-diderot.fr/) 

See also the list of [contributors](https://github.com/samuelmurail/gromacs_py/contributors) who participated in this project.

## License

This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details.
