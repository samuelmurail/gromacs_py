import os
import tools.os_command as os_command

VMD_BIN=os_command.which('vmd_MACOSXX86','vmd','vmd_LINUXAMD64')


#VMD_BIN="/Applications/VMD 1.9.3.app/Contents/vmd/vmd_MACOSXX86"
VMD_MOD_DIRNAME = os.path.dirname(os.path.abspath(__file__))

def insert_mol(pdb_in, pdb_out, out_folder, mol_chain, mol_length, check_file_out = True):
    """
    DEPRECATED, use pdb_manip.insert_mol instead
    
    Use vmd to insert molecules defined by chain ID ``mol_chain`` in a water solvant.

    :param pdb_in: path of input pdb file
    :type pdb_in: str

    :param pdb_out: name of output pdb file
    :type pdb_out: str

    :param out_folder: path of the ouput directory
    :type out_folder: str

    :param mol_chain: chain ID of the molecule to be inserted, 
    :type mol_chain: str

    :param mol_length: number of residue of individual molecule to be inserted
    :type mol_length: int

    :param check_file_out: flag to check or not if file has already been created.
        If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    .. warning::
        The ``pdb_in`` file must contain alredy a concatenated system with a ligand (chain: ``mol_chain``) and a solvated system.
    """

    # Create the out_folder:
    pdb_out = out_folder+"/"+pdb_out
    osCommand.create_dir(out_folder)

    print("\n\nInsert mol in system")

    # Check if output files exist: 
    if check_file_out and os.path.isfile(pdb_out) :
        print("Insert Mol", pdb_out, "already exist")
        return(None)

    script_vmd = VMD_MOD_DIRNAME+"/vmd_script/insert_mol.tcl"
    cmd_insert = os_command.Command([VMD_BIN, "-dispdev", "text", "-e", script_vmd, pdb_in,"-args",mol_chain, str(mol_length), pdb_out])

    cmd_insert.run()

    return(pdb_out)


# import tools.osCommand as osCommand
# osCommand.which('vmd')
