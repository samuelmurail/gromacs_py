#!/usr/bin/env python3
# coding: utf-8

""" Collection of function related to os and sys operations.
"""

import os
import subprocess
import operator

__author__ = "Samuel Murail"

# Test folder path
OS_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(OS_LIB_DIR+"/../test/input/")
TEST_OUT = 'gromacs_py_test_out/os_command/'


def which(*program_list):
    """ find and return the path of a program
    Look for all combination within the `$PATH` env variable

    :param program_list: list of program name
    :type program_list: list of str

    :return: path of the program
    :rtype: str

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
                # Use expanduser in case of ~ caracter
                exe_file = os.path.expanduser(os.path.join(path.replace("\\ ", " "), program))
                #print(exe_file)
                if os.path.isfile(exe_file):
                    return exe_file
    print("program not found !!")
    raise IOError("program not found :", program_list)

def is_exe(fpath):
    """ Check is a file path exist and if user has access to it

    :param fpath: file path
    :type fpath: str

    :return: if the file is an executable
    :rtype: bool

    :Example:

    >>> import tools.os_command as os_command
    >>> print(is_exe('/bin/ls'))
    True

    """

    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def create_and_go_dir(dir_name):
    """ Create the path to a directory and change path in it.

    :param dir_name: directorie name
    :type dir_namet: str

    :Example:

    >>> import tools.os_command as os_command
    >>> start_dir = os.getcwd()
    >>> create_and_go_dir(TEST_OUT+"tmp")
    >>> print("Path: ", os.getcwd()) #doctest: +ELLIPSIS
    Path: .../gromacs_py_test_out/os_command/tmp
    >>> os.chdir(start_dir)
    """

    dir_name = os.path.expanduser(dir_name)
    if not os.path.isdir(dir_name) and dir_name != "":
        os.makedirs(dir_name)
    os.chdir(dir_name)

def create_dir(dir_name):
    """ Create the path to a directory.

    :param dir_name: directorie name
    :type dir_namet: str

    """

    dir_name = os.path.expanduser(dir_name)
    if not os.path.isdir(dir_name) and dir_name != "":
        os.makedirs(dir_name)

def check_file_exist(file):
    """ Check is a file exist.

    :param file: file path
    :type file: str

    :return: if the file exist
    :rtype: bool

    :Example:

    >>> import tools.os_command as os_command
    >>> test_exist = check_file_exist(TEST_PATH+"/1y0m.pdb")
    >>> print("1y0m.pdb exist: ", test_exist)
    1y0m.pdb exist:  True
    """

    file = os.path.expanduser(file)
    return os.path.isfile(file)

def delete_file(file):
    """ Delete a file.

    :param file: file path
    :type file: str

    :return: operation sucess
    :rtype: bool
    """

    file = os.path.expanduser(file)
    return os.remove(file)

def check_file_and_create_path(file):
    """ Check if file exist and create path if not available

    :param file: file path
    :type file: str

    :return: File existance
    :rtype: bool
    """

    file = os.path.expanduser(file)
    if os.path.isfile(file):
        return True
    if os.path.dirname(file) != "":
        create_dir(os.path.dirname(file))
    return False

def full_path_and_check(file):
    """ Return the full path of a file

    :param file: file path
    :type file: str

    :return: File path
    :rtype: str
    """

    if os.path.isfile(file):
        return os.path.abspath(file)

    raise IOError("File cound not be founded :"+file)

def get_directory(file):
    """ Return the path of a file directory

    :param file: file path
    :type file: str

    :return: File path
    :rtype: str
    """

    directory = os.path.dirname(file)
    if directory == "":
        directory = "."
    return directory

class Command:
    """The Command class is a way to launch bash command and mainly gromacs

    :param cmd: command list
    :type cmd: list

    :param env: environment variable
    :type env: dict

    :Example:

    >>> import tools.os_command as os_command
    >>> cmd_list = ['ls','-a',TEST_PATH]
    >>> cmd_test = Command(list_cmd=cmd_list)
    >>> cmd_test.display() #doctest: +ELLIPSIS
    ls -a ...test/input
    >>> return_code = cmd_test.run(out_data=True)
    >>> print(return_code['stdout']) #doctest: +ELLIPSIS
    .
    ..
    1AWR-055_bestene1-mc.pdb
    1dpx.pdb
    1rxz.pdb
    ...
    4n1m.pdb
    ...
    <BLANKLINE>
    """

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
        """ Define the environment of the ``Command`` object.
        """

        self.env = my_env

    def display(self):
        """ Show ``Command`` object that will be launch.
        Show only the name of the command (*eg.* `gmx`) instead of the full path.
        Show relative path for files in the command.
        """

        relative_path_list = []

        for i, arg in enumerate(self.cmd):
            # To avoid showing the full pass of the program extract only the name:
            if i == 0:
                relative_path = os.path.split(arg)[1]
            else:
                try:
                    relative_path = os.path.relpath(arg)
                except:
                    relative_path = arg
            relative_path_list.append(relative_path)
        print(" ".join(relative_path_list))

    def display_raw(self):
        """ Show ``Command`` object that will be launch.
        Show the full path of the command as well as the 
        full path for files in the command.
        """

        print(" ".join(self.cmd))

    def run(self, com_input="", display=False, out_data=False):
        """ Launch ``Command`` object that will be launch.
        return programm output is `out_data` is set to `True`
        """

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
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

    print("-Test os_command module:")

    print("tools.os_command:  \t", doctest.testmod())
