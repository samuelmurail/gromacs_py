# coding: utf-8
import os
import tools.os_command as os_command

#PYMOL_BIN="/Applications/PyMOL.app/Contents/bin/pymol"
#PYMOL_BIN=osCommand.which('pymol')
PYMOL_BIN='pymol'

def make_peptide(sequence, pdb_out, check_file_out = True):
    """
    DEPRECATED, use a pdb_manip method instead.
    Use Pymol to create a linear peptide structure.

    :param sequence: peptide sequence 
    :type sequence: str

    :param pdb_out: name of output pdb file
    :type pdb_out: str

    :param check_file_out: flag to check or not if file has already been created.
        If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    .. warning::
        The ``pdb_in`` file must contain alredy a concatenated system with a ligand (chain: ``mol_chain``) and a solvated system.
    """


    # Create and go in out_folder:
    # This is necessary for the topologie creation
    out_folder = os.path.dirname(pdb_out)
    #print(out_folder)
    os_command.create_dir(out_folder)


    print("-Make peptide using pymol")

    # Check if output files exist: 
    if check_file_out and os_command.check_file_and_create_path(pdb_out) :
    	print("make_peptide", pdb_out, "already exist")
    	return(os.path.abspath(pdb_out))

    pymol_code = "fab G"+sequence+";cmd.alter('resid 1', 'resn=\"ACE\"');cmd.alter('resid 1 and name CA', 'name=\"CH3\"');remove resid 1 and name H+N; sort; save "+pdb_out+", not hydrogens"

    cmd_pep = os_command.Command([PYMOL_BIN, "-d", pymol_code, "-c"])
    #pymol_code = "fab "+sequence+"\n"+"save "+pdb_out+"\n"+"quit \n"

    # The pymol script should be : But the alter command don't work !!!!!
    #fab GGPH
    #remove resid 1 and name H+N
    #alter resid 1, resn="ACE"
    #save test.pdb

    cmd_pep.run()

    pdb_out = os.path.abspath(pdb_out)

    return(pdb_out)

