#!/usr/bin/env python3
# coding: utf-8
#####################################
#########    PDB IN/OUT    ##########
#####################################
__author__ = "Samuel Murail"

import sys
import  os
# Needed for doctest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import tools.osCommand as osCommand


class coor(object):
    """ Topologie base on coordinates like pdb or gro.

    The coor object containt only a dictionnary of atoms indexed on the atom num.


    :param atom_dict: dictionnary of atom
    :type atom_dict: dict

    **Atom dictionnary parameters**
    
    :param field: pdb field
    :type field: str

    :param num: atom number
    :type num: int
 
    :param name: atom name
    :type name: str
 
    :param alter_loc: atom number
    :type alter_loc: str
 
    :param res_name: residue name (3 letters)
    :type res_name: str
 
    :param chain: chain ID
    :type chain: str
 
    :param res_num: residue number (based on pdb file)
    :type res_num: int
 
    :param uniq_resid: unique residue number 
    :type uniq_resid: int
 
    :param insert_res: atom number
    :type insert_res: str
 
    :param x: atom number
    :type x: float
 
    :param y: atom number
    :type y: float
 
    :param z: atom number
    :type z: float
 
    :param occ: occupation
    :type occ: float
 
    :param beta: beta flactor
    :type beta: float


    .. note::
        The atom num index in the dictionnary, is not the same as the ``atom_num`` field of the dictionnary.

    .. note::
        Files necessary for testing : ../test/input/1y0m.pdb, ../test/input/1rxz.pdb and ../test/input/4n1m.pdb. 
        To do the unitary test, execute pdb_mani.py (-v for verbose mode)

    .. todo::
        Add an atom class ?



    """
    def __init__(self):
        self.atom_dict = dict()

    def read_pdb(self, pdb_in, pqr_format = False):
        """Read a pdb file and return atom informations as a dictionnary indexed on the atom num.
        The fonction can also read pqr files if specified with ``pqr_format = True``, it will only change 
        the column format of beta and occ factors.

        :param pdb_in: path of the pdb file to read
        :type pdb_in: str
    
        :param pqr_format: Flag for .pqr file format reading.
        :type pqr_format: bool, default=False
    
        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found

        """  

        atom_index = 0
        uniq_resid = -1
        old_res_num = -1
    
        with open(pdb_in) as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' or line[:6] == "HETATM":
    
                    field       = line[:6].strip()
                    atom_num    = int(line[6:11])
                    atom_name   = line[12:16].strip()
                    alter_loc   = line[16:17]
                    res_name    = line[17:20].strip()
                    chain       = line[21]
                    res_num     = int(line[22:26])
                    insert_res  = line[26:27]
                    x, y, z     = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    if pqr_format:
                        occ, beta   = line[54:62].strip(), line[62:70].strip()
                    else:
                        occ, beta   = line[54:60].strip(), line[60:66].strip()
    
                    if occ == "":
                        occ = 0.0
                    else:
                        occ = float(occ)
    
                    if beta == "":
                        beta = 0.0
                    else:
                        beta = float(beta)
    
                    if res_num != old_res_num:
                        uniq_resid +=1
                        old_res_num = res_num
    
                    atom = {"field" : field, 
                        "num"       : atom_num, 
                        "name"      : atom_name, 
                        "alter_loc" : alter_loc, 
                        "res_name"  : res_name, 
                        "chain"     : chain, 
                        "res_num"   : res_num,
                        "uniq_resid": uniq_resid,
                        "insert_res": insert_res,
                        "x"         : x,
                        "y"         : y,
                        "z"         : z,
                        "occ"       : occ,
                        "beta"      : beta}
    
                    self.atom_dict[atom_index] = atom
                    atom_index +=1
        print("Succeed to read file",pdb_in,", ",atom_index,"atoms found")

    def write_pdb(self, pdb_out):
        """Write a pdb file.

        :param pdb_out: path of the pdb file to write
        :type pdb_out: str

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.write_pdb('../test/output/pdb_manip_test/tmp.pdb')
        Succeed to save file ../test/output/pdb_manip_test/tmp.pdb

        """  
    
        filout = open(pdb_out, 'w')
        
        for atom_num, atom in sorted(self.atom_dict.items()):
            #print(pdb_dict[atom_num]["name"])
            filout.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                atom["field"],
                atom["num"], 
                atom["name"], 
                atom["alter_loc"], 
                atom["res_name"], 
                atom["chain"], 
                atom["res_num"], 
                atom["insert_res"], 
                atom["x"], 
                atom["y"], 
                atom["z"], 
                atom["occ"], 
                atom["beta"]))
        filout.write("TER\n")
        filout.close()

        print("Succeed to save file",pdb_out)


    def get_aa_seq(self):
        """Get the amino acid sequence from a coor object.
        
        :return: dictionnary of sequence indexed by the chain ID
        :rtype: dict

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSVVKVLFDYKVQREDELTFTKSVIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        .. warning::
            If atom chains are not arranged sequentialy (A,A,A,B,B,A,A,A ...), the first atom seq will be overwritten by the last one. 

        """  
        
        aa_dict = { "GLY":"G", 
                    "HIS":"H", 
                    "HSE":"H", 
                    "HSD":"H", 
                    "HSP":"H",
                    "ARG":"R",
                    "LYS":"K",
                    "ASP":"D",
                    "GLU":"E",
                    "SER":"S",
                    "THR":"T",
                    "ASN":"N",
                    "GLN":"Q",
                    "CYS":"C",
                    "SEC":"U",                   
                    "PRO":"P",
                    "ALA":"V",                        
                    "ILE":"I",
                    "PHE":"F",
                    "TYR":"Y",
                    "TRP":"W",
                    "VAL":"V",
                    "LEU":"L",
                    "MET":"M"
                    }

                                                
        # Get CA atoms 
        CA_index_list = self.get_index_selection( {"name" : ["CA"]} )

        seq = ""
        seq_dict = {}
        chain_first = self.atom_dict[CA_index_list[0]]['chain']
        
        for index in sorted(CA_index_list):
            loop_atom = self.atom_dict[index]
            
            if loop_atom['chain'] != chain_first:
                seq_dict[chain_first] = seq
                seq = ""
                chain_first = loop_atom['chain']
                
            seq = seq + aa_dict[loop_atom['res_name']]
            
        seq_dict[chain_first] = seq

        return(seq_dict)

    def get_aa_num(self):
        """Get the amino acid number a coor object.
        
        :return: Number of residues
        :rtype: int

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61

        .. note::
            Only count Ca atoms, this may not be the best choice ?    
        
        """  

        CA_index_list = self.get_index_selection( {"name" : ["CA"]} )

        return( len(CA_index_list) )


    def change_pdb_field(self, change_dict):
        """Change all atom field of a coor object, the change is based on the change_dict dictionnary.
        
        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict 

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSVVKVLFDYKVQREDELTFTKSVIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.change_pdb_field(change_dict = {"chain" : "B"}) #doctest: +ELLIPSIS
        <tools.pdb_manip.coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'B': 'TFKSVVKVLFDYKVQREDELTFTKSVIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        """  
    
        for atom_num, atom in self.atom_dict.items():
            for change, val in change_dict.items():
                atom[change] = val
    
        return(self)

    def change_index_pdb_field(self, index_list, change_dict):
        """Change all atom field of a part of coor object defined by ``index``, the change is based on the change_dict dictionnary.
        
        :param index_list: list of atom index to change 
        :type index_list: list  

        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict 

        :Example:
        
        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> res_826_852 = prot_coor.get_index_selection({'res_num' : range(826,852)})
        >>> prot_coor.change_index_pdb_field(index_list = res_826_852, change_dict = {"chain" : "B"}) #doctest: +ELLIPSIS
        <tools.pdb_manip.coor object at ...>
        >>> prot_seq = prot_coor.get_aa_seq()
        >>> prot_seq == {'A': 'TFKSVVKVLFDYKVQREDELTFTKSVIIQNVEKQD', 'B': 'GGWWRGDYGGKKQLWFPSNYVEEMIN'}
        True

        """  

        for atom_num in index_list:
            for change, val in change_dict.items():
                self.atom_dict[atom_num][change] = val
    
        return(self)
    

    def select_part_dict(self, selec_dict):
        """Select atom of a coor object defined, the selection is based on the change_dict dictionnary.
        Return a new coor object.
        
        :param selec_dict: change ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict 

        :return: a new coor object
        :rtype: coor

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61
        >>> prot_20_coor = prot_coor.select_part_dict(selec_dict = {'res_num' : list(range(791,800))})
        >>> prot_20_coor.get_aa_seq()
        {'A': 'TFKSVVKVL'}
        >>> prot_20_coor.get_aa_num()
        9

        """  

        coor_out = coor()
    
        for atom_num, atom in self.atom_dict.items():
            selected = True
            for selection in selec_dict.keys():
                #print("select",selection, selec_dict[selection],".")
                #print("atom",atom[selection],".")
                if  atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                coor_out.atom_dict[atom_num] = atom
    
        return(coor_out)
    
    def get_index_selection(self, selec_dict):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of index of selected atoms.
        
        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict 

        :return: list of atom index
        :rtype: list of int

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_index_selection({'res_num' : [826,827]})
        [297, 298, 299, 300, 301, 302, 303, 304]

        """  
    
        index_list = []
    
        for atom_num, atom in self.atom_dict.items():
            selected = True
            #print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                #print("select",selection, selec_dict[selection],".")
                #print("selection=",selection)
                #print("atom:",atom)
                #print("atom",atom[selection],".")
                if  atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                #print(atom)
                index_list.append(atom_num)
    
        return( index_list )
    
    def get_uniq_res_selection(self, selec_dict):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of unique residue of selected atoms.
        
        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict 

        :return: list of unique residue
        :rtype: list of int

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_uniq_res_selection({'res_num' : [826,827]})
        [35, 36]

        .. note::
            Unique residues are assigned sequentialy (0, 1, 2 ...) when reading the pdb file. As res_num field is extracted from the pdb file field.

        """  

        uniq_res_list = []
    
        for atom_num, atom in self.atom_dict.items():
            selected = True
            for selection in selec_dict.keys():
                #print("select",selection, selec_dict[selection],".")
                #print("atom",atom[selection],".")
                if  atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                #print(atom)
                uniq_res_list.append(self.atom_dict[atom_num]["uniq_resid"])
    
        return( list(set(uniq_res_list)) )

    def del_atom_index(self, index_list):
        """Delete atoms of a coor object defined by their ``index``.
        
        :param index_list: list of atom index to delete
        :type index_list: list  

        :Example:
        
        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSVVKVLFDYKVQREDELTFTKSVIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> res_810_852 = prot_coor.get_index_selection({'res_num' : range(810,852)})
        >>> prot_coor.del_atom_index(index_list = res_810_852) #doctest: +ELLIPSIS
        <tools.pdb_manip.coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSVVKVLFDYKVQREDE'}

        """  

        for index in index_list:
            del self.atom_dict[index]

        return(self)
    
    def correct_chain(self, Ca_cutoff = 4.5):
        """Correct the chain ID's of a coor object, by checking consecutive Calphas atoms distance.

        
        :param Ca_cutoff: cutoff for distances between Calphas atoms (Å)
        :type Ca_cutoff: float, default=4.5

        :Example:
        
        >>> import tools.pdb_manip as pdb_manip
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/input/1y0m.pdb')
        Succeed to read file ../test/input/1y0m.pdb ,  648 atoms found
        >>> res_810 = prot_coor.get_index_selection({'res_num' : [810]})
        >>> prot_coor = prot_coor.del_atom_index(index_list = res_810)
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSVVKVLFDYKVQREDETFTKSVIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.correct_chain() #doctest: +ELLIPSIS
        Chain: A  Residue: 0 to 18
        Chain: B  Residue: 20 to 60
        <tools.pdb_manip.coor object at ...>
        >>> # As a residue is missing, Calphas after residue 18 is no more consecutive


        .. note::
            This is specially usefull for pdb2gmx which cut the protein chains based on the chain ID's

        """  

        """ Check distance between consecutive Caphas
        If they are away of Ca_cutoff, they are considered 
        as in different chain"""
    
        Ca_atom = self.select_part_dict({"name" : ["CA"]})
        first_flag = True
        chain_res_list = []
        res_list = []
    
        # Identify Chain uniq_resid
        # Need to use sorted to be sure to check consecutive residues (atoms)
        for key, atom in sorted(Ca_atom.atom_dict.items()):
            if first_flag:
                first_flag = False
            else:
                distance = coor.atom_dist(atom, old_atom)
                #print(distance)
                if distance < Ca_cutoff :
                    old_atom = atom
                else:
                    #print("New chain")
                    chain_res_list.append(res_list)
                    res_list = []
            res_list.append(atom['uniq_resid'])
            old_atom = atom
        chain_res_list.append(res_list)
        
        # Change chain ID :

        #print(Ca_atom.atom_dict)

        for i, chain_res in  enumerate(chain_res_list):
            print("Chain:",chr(65+i)," Residue:",chain_res[0],"to",chain_res[-1])
            chain_index = self.get_index_selection({'uniq_resid' : chain_res})
            #print(chain_index)
            self.change_index_pdb_field(chain_index, {"chain" : chr(65+i)} )
            #pdb_dict_out.update(chain_dict)
    
        #print(pdb_dict_out)
    
        return(self)
    
    def correct_HIS_name(self):
        """ Get his protonation state from pdb2pqr and replace HIS resname 
        with HSE, HSD, HSP resname.
        To do after pdb2pqr, in order that protonation is recognize by pdb2gmx  

        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> import tools.pdb2pqr as pdb2pqr
        >>> # Compute protonation with pdb2pqr: 
        >>> pdb2pqr.compute_pdb2pqr('../test/input/4n1m.pdb', '../test/output/pdb_manip_test/4n1m.pqr')
        Succeed to read file ../test/input/4n1m.pdb ,  2530 atoms found
        Succeed to save file ../test/output/pdb_manip_test/tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain ../test/output/pdb_manip_test/tmp_pdb2pqr.pdb ../test/output/pdb_manip_test/4n1m.pqr
        0
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/output/pdb_manip_test/4n1m.pqr', pqr_format = True)
        Succeed to read file ../test/output/pdb_manip_test/4n1m.pqr ,  2548 atoms found
        >>> HSD_index = prot_coor.get_index_selection({'res_name' : ['HSD'], 'name':['CA']})
        >>> print(len(HSD_index))
        5
        >>> HSE_index = prot_coor.get_index_selection({'res_name' : ['HSE'], 'name':['CA']})
        >>> print(len(HSE_index))
        0
        >>> HSP_index = prot_coor.get_index_selection({'res_name' : ['HSP'], 'name':['CA']})
        >>> print(len(HSP_index))
        0
        >>> prot_coor.correct_HIS_name() #doctest: +ELLIPSIS
        <tools.pdb_manip.coor object at 0x...
        >>> HIS_index = prot_coor.get_index_selection({'res_name' : ['HIS'], 'name':['CA']})
        >>> print(len(HIS_index))
        0


        .. note::
            This function seems useless. Since last version of pdb2pqr residue name seems correct.

        """  
    
        # FIND HISTIDINE res 
    
        # HSD:
        hsd_uniq_res = self.get_uniq_res_selection({"res_name" : ["HIS"], "name" : ["HD1"]})
        # HSE:
        hse_uniq_res = self.get_uniq_res_selection({"res_name" : ["HIS"], "name" : ["HE2"]})
        # HSP: find res in common with both hsd and hse
        hsp_uniq_res = [res for res in hsd_uniq_res if res in hse_uniq_res]
        # remove HSP res from HSE HSD list
        if len(hsp_uniq_res) != 0:
            for res in hsp_uniq_res:
                hsd_uniq_res.remove(res)
                hse_uniq_res.remove(res)
    
        # Replace HIS resname :
        all_his_uniq_res = hsd_uniq_res + hse_uniq_res + hsp_uniq_res
    
        for atom_num, atom in self.atom_dict.items():
            selected = True
            if atom["uniq_resid"] in all_his_uniq_res:
                if atom["uniq_resid"] in hsd_uniq_res:
                    atom["res_name"] = "HSD"
                elif atom["uniq_resid"] in hse_uniq_res:
                    atom["res_name"] = "HSE"
                else :
                    atom["res_name"] = "HSP"
                #print(atom)
        
        return(self)
    
    def correct_CYS_name(self):
        """ Correct the CYS resname from pdb2pqr
    
        :Example:

        >>> import tools.pdb_manip as pdb_manip
        >>> import tools.pdb2pqr as pdb2pqr
        >>> # Compute protonation with pdb2pqr: 
        >>> pdb2pqr.compute_pdb2pqr('../test/input/1dpx.pdb', '../test/output/pdb_manip_test/1dpx.pqr')
        Succeed to read file ../test/input/1dpx.pdb ,  1192 atoms found
        Succeed to save file ../test/output/pdb_manip_test/tmp_pdb2pqr.pdb
        pdb2pqr.py --ff CHARMM --ffout CHARMM --chain ../test/output/pdb_manip_test/tmp_pdb2pqr.pdb ../test/output/pdb_manip_test/1dpx.pqr
        0
        >>> prot_coor = pdb_manip.coor()
        >>> prot_coor.read_pdb('../test/output/pdb_manip_test/1dpx.pqr', pqr_format = True)
        Succeed to read file ../test/output/pdb_manip_test/1dpx.pqr ,  1960 atoms found
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['ISU']})
        >>> print(len(Isu_index))
        16
        >>> prot_coor.correct_CYS_name() #doctest: +ELLIPSIS
        <tools.pdb_manip.coor object at 0x...
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['ISU']})
        >>> print(len(Isu_index))
        0

        """

        # FIND ISU res 
        ISU_index_list = self.get_index_selection({"res_name" : ["ISU"]})
        if len(ISU_index_list) == 0:
            #print("Nothing to Fix")
            return(self)
    
        # Replace CYS resname :
    
        for atom_num in ISU_index_list:
            selected = True
            self.atom_dict[atom_num]["res_name"] = "CYS"
            if self.atom_dict[atom_num]["name"] == "1CB":
                self.atom_dict[atom_num]["name"] = "CB"
            if self.atom_dict[atom_num]["name"] == "1SG":
                self.atom_dict[atom_num]["name"] = "SG"
    
        return(self)
      
    
    def insert_mol(self, pdb_out, out_folder, mol_chain, mol_length, check_file_out = True):
        """
        Insert molecules defined by chain ID ``mol_chain`` in a water solvant.
    
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
    
        # Select prot atoms :
        prot_CA = self.select_part_dict(selec_dict = {'name' : ['CA']})  
        water_O = self.select_part_dict(selec_dict = {'res_name' : ['SOL'], 'name':['OW']})  
        
        print(len(prot_CA.atom_dict), len(water_O.atom_dict)) 

        return(pdb_out)


    @staticmethod
    def atom_dist(atom_1, atom_2):
        """Compute the distance between 2 atoms.

        :param atom_1: atom dictionnary
        :type atom_1: dict

        :param atom_2: atom dictionnary
        :type atom_2: dict

        :Example:
        
        >>> import tools.pdb_manip as pdb_manip
        >>> atom_1 = {'x':0.0, 'y':0.0, 'z': 0.0}
        >>> atom_2 = {'x':0.0, 'y':1.0, 'z': 0.0}
        >>> atom_3 = {'x':1.0, 'y':1.0, 'z': 1.0}
        >>> coor.atom_dist(atom_1, atom_2)
        1.0
        >>> coor.atom_dist(atom_1, atom_3)
        1.7320508075688772

        """  


        distance = ( (atom_1['x']-atom_2['x'])**2 + (atom_1['y']-atom_2['y'])**2 + (atom_1['z']-atom_2['z'])**2 ) ** 0.5
        return(distance)
    
    
    @staticmethod
    def concat_pdb(*pdb_in_files, pdb_out):
        """Concat a list of pdb files in one.

        :param pdb_in_files: list of pdb files
        :type pdb_in_files: list

        :param pdb_out: atom dictionnary
        :type pdb_out: dict

        :Example:
        
        >>> import tools.pdb_manip as pdb_manip
        >>> coor.concat_pdb('../test/input/1y0m.pdb','../test/input/1rxz.pdb', pdb_out = '../test/output/pdb_manip_test/tmp_2.pdb')
        Concat : ../test/input/1y0m.pdb
        Concat : ../test/input/1rxz.pdb
        Succeed to save ../test/output/pdb_manip_test/tmp_2.pdb

        """  


        if osCommand.check_file_and_create_path(pdb_out):
            print("File "+pdb_out+" already exist")
            return
    
        filout = open(pdb_out, 'w')
        count = 0
    
        for pdb_in in pdb_in_files:
            print("Concat :",pdb_in)
            with open(pdb_in) as pdbfile:
                for line in pdbfile:
                    if (count == 0 and line[:6] ==  "CRYST1") or line[:4] == 'ATOM' or line[:6] == "HETATM":
                        filout.write(line)
            count += 1
    
        filout.close()
        print("Succeed to save",pdb_out)


if __name__ == "__main__":

    import doctest
    import shutil
    shutil.rmtree('../test/output/pdb_manip_test', ignore_errors=True)
    doctest.testmod("pdb_manip")
    # Erase all test files
    shutil.rmtree('../test/output/pdb_manip_test', ignore_errors=True)

