#!/usr/bin/env python3
# coding: utf-8

#####################################
#########     SHORTCUTS    ##########
#####################################

import os
import subprocess
import operator



def which(*program_list):
    """ find and return the path of a programm
    Look for all combination with the PATH env variable"""

    for program in program_list:
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                #print(path)
                # Use expanduser in case of ~ caracter
                exe_file = os.path.expanduser(os.path.join(path.replace("\ "," "), program))
                #print(exe_file)
                if os.path.isfile(exe_file):
                    return exe_file
    print("program not found !!")
    raise IOError("program not found :",program_list)

def create_or_go_dir(dir_name):
    dir_name = os.path.expanduser(dir_name)
    if not os.path.isdir(dir_name) and dir_name != "":
        os.makedirs(dir_name)
    os.chdir(dir_name) 

def create_dir(dir_name):
    dir_name = os.path.expanduser(dir_name)
    if not os.path.isdir(dir_name) and dir_name != "":
        os.makedirs(dir_name)

def check_file_exist(file):
    file = os.path.expanduser(file)
    if os.path.isfile(file):
        return True

def delete_file(file):
    file = os.path.expanduser(file)
    if os.remove(file):
        return True

def check_file_and_create_path(file):
    """ Check if file exist and create dir if not available """
    file = os.path.expanduser(file)
    if os.path.isfile(file):
        return True
    elif os.path.dirname(file) != "":
        create_dir(os.path.dirname(file))
    return False

def full_path_and_check(file):
    if os.path.isfile(file):
        return os.path.abspath(file)
    else:
        raise IOError("File cound not be founded :"+file)

def get_directory(file):
    directory = os.path.dirname(file)
    if directory == "":
        directory = "."
    return(directory)

class command(object):

    def __init__(self, list_cmd, my_env = None, **kwargs):
        self.cmd = list_cmd
        self.env = my_env
        # Add supplementary argument for the command
        # Mainly usefull for gromacs
        # Add "-" to the key, eg ter -> -ter
        if kwargs is not None:
            # Use sorted to have same order of command in doctest
            for key, value in sorted( kwargs.items(), key=operator.itemgetter(0)):
                #print(key, value )
                self.cmd.append("-"+key)
                self.cmd.append(value)
        #print("Cmd:",self.cmd)

        #sorted(x.items(), key=operator.itemgetter(1))

    def define_env(self, my_env):
        self.env = my_env   

    def display(self):
        
        relative_path_list = []
        
        for arg in self.cmd:
            #print(arg)
            try: 
                relative_path = os.path.relpath(arg)
            except:
                relative_path = arg
            relative_path_list.append(relative_path)
        print(" ".join(relative_path_list))

    def display_raw(self):
        print(" ".join(self.cmd))

    def run(self, com_input = "", display = False, out_data = False):

        proc = subprocess.Popen(self.cmd, 
            stdin=subprocess.PIPE, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            env = self.env)
        (stdout_data, stderr_data) = proc.communicate(com_input.encode())

        if display:
            print(self.display())
            print(stdout_data.decode('utf-8'))
            print(stderr_data.decode('utf-8'))


        # Check if command is successfull
        if proc.returncode == 0:
            if out_data:
                return( {'returncode':proc.returncode, 'stdout':stdout_data.decode('utf-8') , 'stderr':stderr_data.decode('utf-8') } )
            else:
                return( proc.returncode )
        else:
            print("The following command could not be executed correctly :")
            self.display()
            print(stdout_data.decode('utf-8'))
            print(stderr_data.decode('utf-8'))
            raise Error()

    def slurm_peptide(self, script_name, pep_list, output="slumr_pep.out"):
        """
        #!/bin/bash
        #SBATCH --output=test.out
        #SBATCH --array=1-5
        #SBATCH -p production
        case "$SLURM_ARRAY_TASK_ID" in 
           1) fichier='/scratch/user/rey/protein_1';;
           2) fichier='/scratch/user/rey/protein_2';;
           3) fichier='/scratch/user/rey/protein_3';;
           4) fichier='/scratch/user/rey/protein_4';;
           5) fichier='/scratch/user/rey/protein_5';;
        esac
        drun opendocking babel -i ${fichier}.pdb -o ${fichier}.mol2
        """

        filout = open(script_name, 'w')
        
        filout.write("#!/bin/bash\n")
        filout.write("#SBATCH --output={}\n".format(output) )
        filout.write("#SBATCH --array=1-{}\n".format(len(pep_list)))
        filout.write("#SBATCH -p production\n")
        filout.write("case \"$SLURM_ARRAY_TASK_ID\" in\n")
        for i, pep in enumerate(pep_list):
            filout.write("{}) peptide=\'{}\';;\n".format(i,pep))
        filout.write("esac\n")
        filout.write("drun gromacs_lib sim_pep_prot.py -seq ${peptide}.pdb -npep 13 -Pep_time 1 ")
        filout.write("-fsys 1y0m_water_ions/03_prod/prod_1y0m.gro -psys 1y0m_water_ions/top/1y0m_water_ions_water_ion.top -ssys 1y0m_water_ions/03_prod/prod_1y0m.tpr \n")
        filout.write("-o ${peptide} -n ${peptide} -dt 0.005 -em_steps 5000 -dt_HA 0.002 -HA_time 0.25 -CA_time 1 -CA_LOW_time 5 -PROD_time 500 \n")
        filout.close()


    ## General args:
    #parser.add_argument('-o', action="store", dest="o", help='Output Directory', required=True)
    #parser.add_argument('-n', action="store", dest="name", help='Output file name', required=True)
    #parser.add_argument('-dt', action="store", dest="dt", help='Equilibration dt, default=0.005 (5 fs)', type=float, default=0.005)
    #parser.add_argument('-em_steps', action="store", dest="em_steps", help='Minimisation steps, default = 5000', type=int, default=5000)
    ## Equilibration and production args:
    #parser.add_argument('-dt_HA', action="store", dest="dt_HA", help='Equilibration dt, default=0.005 (5 fs)', type=float, default=0.002)
    #parser.add_argument('-HA_time', action="store", dest="HA_time", help='Equilibration with HA constraint time(ns), default = 0.25ns', type=float, default=0.25)
    #parser.add_argument('-CA_time', action="store", dest="CA_time", help='Equilibration with HA constraint time(ns), default = 1ns', type=float, default=1)
    #parser.add_argument('-CA_LOW_time', action="store", dest="CA_LOW_time", help='Equilibration with HA constraint time(ns), default = 5ns', type=float, default=5)
    #parser.add_argument('-PROD_time', action="store", dest="Prod_time", help='Production time(ns), default = 100ns', type=float, default=100)
    ## mdrun args:
    #parser.add_argument('-nt', action="store", dest="nt", help='Total number of threads to start, default=0', type=float, default=0)
    #parser.add_argument('-ntmpi', action="store", dest="ntmpi", help='Number of thread-MPI threads to start, default=0', type=float, default=0)
    #parser.add_argument('-gpu_id', action="store", dest="gpuid", help='List of GPU device id-s to use, default=\"\" ', default="None")

#sim_pep_prot.py -seq ARG -npep 10 -Pep_time 0.1 -fsys 1y0m_water_ions/03_prod/prod_1y0m.gro -psys 1y0m_water_ions/top/1y0m_water_ions_water_ion.top -o test_prot_AR -n prot_AR




if __name__ == '__main__':

    test = command(["ls","-lsh"])
    print(test)
    test.display()
    out = test.run_redirect_outerr()
    print(out)
    
    test.run()