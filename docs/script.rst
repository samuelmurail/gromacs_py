Script
=======================================

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
	