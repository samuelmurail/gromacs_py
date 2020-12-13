Script
=======================================

***************************************
Tutorial
***************************************

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
   create_top.py -f gromacs_py/test_files/1y0m.pdb -o tmp/1y0m/top -vsite

   # Minimize the protein structure
   minimize_pdb.py -f tmp/1y0m/top/1y0m_pdb2gmx_box.pdb -p tmp/1y0m/top/1y0m_pdb2gmx.top -o tmp/1y0m/em/  -n em_1y0m -nt 2

   # Add water and ions
   solvate_ions.py -f tmp/1y0m/em/em_1y0m_compact.pdb -p tmp/1y0m/top/1y0m_pdb2gmx.top -o tmp/1y0m_water_ions/top/  -n 1y0m_water_ions

   # Minimize the system
   minimize_pdb.py -f tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.gro -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/em/  -n em_1y0m

   # Do three small equilibrations with postion contraints on heavy atoms (first), Carbon alpha (second) and low constraint on Carbon alpha (third)
   equi_3_step.py -f tmp/1y0m_water_ions/em/em_1y0m_compact.pdb -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/  -n 1y0m -HA_time 0.1 -dt_HA 0.002 -CA_time 0.1 -CA_LOW_time 0.1 -dt 0.004 -maxwarn 1

   # Small production run of 0.1 ns
   production.py -f tmp/1y0m_water_ions/02_equi_CA_LOW/equi_CA_LOW_1y0m.gro -p tmp/1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o tmp/1y0m_water_ions/03_prod -n 1y0m -time 0.1 -dt 0.004 -maxwarn 1

   # Extension of the simulation
   extend.py -s tmp/1y0m_water_ions/03_prod/prod_1y0m.tpr -time 0.2

   # Remove simulation files
   rm -r ./tmp

Or simply use one command to do all previous commands:

.. code-block:: bash

   top_em_equi_3_step_prod.py -f gromacs_py/test/input/1y0m.pdb -o tmp/1y0m -vsite -HA_time 0.1 -CA_time 0.1 -CA_LOW_time 0.1 -prod_time 0.3

***************************************
Topologie related:
***************************************

Create topologie
---------------------------------------

.. argparse::
   :filename: ../bin/create_top.py
   :func: parser_input
   :prog: create_top.py
   :nodefault:


Solvate a system
---------------------------------------

.. argparse::
   :filename: ../bin/solvate_ions.py
   :func: parser_input
   :prog: solvate_ions.py
   :nodefault:

***************************************
Simulation:
***************************************

Energy minimization
---------------------------------------

.. argparse::
   :filename: ../bin/minimize_pdb.py
   :func: parser_input
   :prog: minimize_pdb.py
   :nodefault:

Equilibration
---------------------------------------

.. argparse::
   :filename: ../bin/equi_3_step.py
   :func: parser_input
   :prog: equi_3_step.py
   :nodefault:

Production
---------------------------------------

.. argparse::
   :filename: ../bin/production.py
   :func: parser_input
   :prog: production.py
   :nodefault:


***************************************
Peptide related:
***************************************

Create peptide
---------------------------------------

.. argparse::
   :filename: ../bin/create_peptide.py
   :func: parser_input
   :prog: create_peptide.py
   :nodefault:


Minimize cyclic peptide
---------------------------------------

.. argparse::
   :filename: ../bin/minimize_pdb_and_cyclic.py
   :func: parser_input
   :prog: minimize_pdb_and_cyclic.py
   :nodefault:

Insert n copy of a peptide in a system
---------------------------------------

.. argparse::
   :filename: ../bin/sim_pep_prot.py
   :func: parser_input
   :prog: sim_pep_prot.py
   :nodefault:

..
	
	~/Documents/repository/gromacs_py/equi_3_step.py -f 1y0m_water_ions/em/em_1y0m_compact.pdb -p 1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o 1y0m_water_ions/	equi/  -n 1y0m -HA_time 0.1 -CA_time 0.2 -CA_LOW_time 0.4
	