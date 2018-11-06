#!/usr/bin/env python3
# coding: utf-8

#####################################
#########     SHORTCUTS    ##########
#####################################

import os
import subprocess
import operator



def which(*program_list):
    """ find and return the path of a program
    Look for all combination within the `$PATH` env variable

    :param program_list: list of program name
    :type program_list: list of str

    :Example:

    >>> import tools.os_command as os_command
    >>> ls_path = os_command.which('ls')
    >>> print(ls_path)
    /bin/ls

    """

    for program in program_list:
        fpath, fname = os.path.split(program)
        #print(fpath, fname)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                #print(path)
                # Use expanduser in case of ~ caracter
                exe_file = os.path.expanduser(os.path.join(path.replace("\\ ", " "), program))
                #print(exe_file)
                if os.path.isfile(exe_file):
                    return exe_file
    print("program not found !!")
    raise IOError("program not found :", program_list)

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

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
    return os.path.isfile(file)

def delete_file(file):
    file = os.path.expanduser(file)
    return os.remove(file)

def check_file_and_create_path(file):
    """ Check if file exist and create dir if not available 
    """

    file = os.path.expanduser(file)
    if os.path.isfile(file):
        return True
    if os.path.dirname(file) != "":
        create_dir(os.path.dirname(file))
    return False

def full_path_and_check(file):
    
    if os.path.isfile(file):
        return os.path.abspath(file)
    
    raise IOError("File cound not be founded :"+file)

def get_directory(file):
    directory = os.path.dirname(file)
    if directory == "":
        directory = "."
    return directory

class Command:

    def __init__(self, list_cmd, my_env=None, **kwargs):
        self.cmd = list_cmd
        self.env = my_env
        # Add supplementary argument for the command
        # Mainly usefull for gromacs
        # Add "-" to the key, eg ter -> -ter
        if kwargs is not None:
            # Use sorted to have same order of command in doctest
            for key, value in sorted(kwargs.items(), key=operator.itemgetter(0)):
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

    def run(self, com_input="", display=False, out_data=False):

        proc = subprocess.Popen(self.cmd,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                env=self.env)

        (stdout_data, stderr_data) = proc.communicate(com_input.encode())

        if display:
            print(self.display())
            print(stdout_data.decode('utf-8'))
            print(stderr_data.decode('utf-8'))

        # Check if command is successfull
        if proc.returncode == 0:
            if out_data:
                return {'returncode':proc.returncode,
                        'stdout':stdout_data.decode('utf-8'),
                        'stderr':stderr_data.decode('utf-8')}
            return proc.returncode

        print("The following command could not be executed correctly :")
        self.display()
        print(stdout_data.decode('utf-8'))
        print(stderr_data.decode('utf-8'))
        raise Error()


if __name__ == "__main__":

    import doctest
    #import shutil

    print("-Test os_command module:")

    doctest.testmod()

