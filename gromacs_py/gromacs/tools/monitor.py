#!/usr/bin/env python3

# coding: utf-8

""" Collection of function to monitor a simulation in real time.
"""

__author__ = "Samuel Murail"

import matplotlib
import matplotlib.pyplot as plt
import time
import os
from . import os_command

# Test folder path
MONITOR_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(MONITOR_LIB_DIR, "../../test/input/"))


def simulation_plot(proc, func_input_dict, refresh_time=1.0):
    """ This function is used for monitoring a simulation in real time.
    Function can be excecuted by the gromacs.tools.os_command.run_background() function.
    The function monitors a trajectory file, and launch the analysis if the file has been modified.
    It can plot as function of time an analysis of a simulation.
    Analysis is passed as input function.
    """

    file_to_check = func_input_dict[func_input_dict['file_check_ext']]
    function_list = func_input_dict['extract_func']
    time_modif = None
    file_time = None
    count = 1
    x_list = []
    y_list = []

    ###################
    # Set up the plot #
    ###################

    # Remove the matplotlib window buttons:
    matplotlib.rcParams['toolbar'] = 'None'
    # Create the plot in interactive mode
    plt.ion()
    # fig = plt.figure()
    fig, axarr = plt.subplots(len(function_list), sharex=True)

    for i, function in enumerate(function_list):
        y_list.append([])
        x_list.append([])
        axarr[i].set_xlabel('time (ns)')
        axarr[i].set_ylabel(function['term'])
        axarr[i].plot(x_list[i], y_list[i],
                      'ko-', markersize=2,
                      linewidth=0.5,
                      color='blue')
    # show the window
    # figure will be in foreground, but the user may move it to background
    fig.show()
    fig.canvas.set_window_title(file_to_check[:-4])

    while proc.poll() is None:

        time.sleep(refresh_time)
        count += 1

        if os_command.check_file_exist(file_to_check):
            file_time = os.stat(file_to_check).st_mtime

        if time_modif != file_time:

            time_modif = file_time

            for i, function in enumerate(function_list):
                anal = function['func'](func_input_dict)
                # print(anal)

                try:
                    x_list[i].append(anal['time'])
                    y_list[i].append(anal[function['term']])

                    axarr[i].lines[0].set_data(x_list[i], y_list[i])  # set plot data
                    axarr[i].relim()                  # recompute the data limits
                    axarr[i].autoscale_view()         # automatic axis scaling
                except KeyError:
                    print('Energy could not be extract, simulation is probably finished.')

            fig.canvas.flush_events()   # update the plot and take care of window events (like resizing etc.)


def extract_log_dict(func_input_dict, tail_line_num=20):
    """ Read last line of a gromacs ``.log`` file and return a dictionnary
    containing ``time``,  ``step`` and all log energetic terms.
    """

    log_to_check = func_input_dict['log']
    tail_text = os.popen('tail -n {} {}'.format(tail_line_num, log_to_check)).read()

    split_text = tail_text.split('\n')
    log_dict = {}
    ener_read = False
    ener_done = False
    field_len = 15

    i = 0
    while i < len(split_text) and not ener_done:
        line = split_text[i]
        # Find Step Time line:
        if line.strip().startswith('Step'):
            line_split = split_text[i + 1].split()
            log_dict['step'] = int(line_split[0])
            log_dict['time'] = float(line_split[1])
            # print(log_dict)
            # Skip next line (already extracted with time and step)
            i += 2
            continue
        if line.strip().startswith('Energies'):
            ener_read = True
        elif ener_read and len(line) == 0:
            ener_done = True
        elif ener_read:
            next_line = split_text[i + 1]
            for j in range(5):
                field = line[field_len * j:field_len * (j + 1)].strip().replace(" ", "_")
                if len(field) > 0:
                    value = float(next_line[field_len * j:field_len * (j + 1)].strip())
                    log_dict[field] = value
            i += 2
            continue
        i += 1

    return(log_dict)


def print_log_file(proc, func_input_dict, refresh_time=0.01, tail_line_num=20):
    """ Monitor ``.log`` file information.

    Example:

    >>> TEST_OUT = str(getfixture('tmpdir'))
    >>> import sys
    >>> #print(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))
    >>> sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    >>> import gromacs.gmx5 as gmx #doctest: +ELLIPSIS
    Gromacs version is ...
    FORCEFIELD_PATH = ...
    >>> prot = gmx.GmxSys(name='1y0m', coor_file=TEST_PATH+'/1y0m.pdb')
    >>> ###################################
    >>> ####   Create the topologie:   ###
    >>> ###################################
    >>> prot.prepare_top(out_folder=os.path.join(TEST_OUT, 'top_SH3')) #doctest: +ELLIPSIS
    Succeed to read file .../test/input/1y0m.pdb ,  648 atoms found
    Succeed to save file tmp_pdb2pqr.pdb
    pdb2pqr.py --ff CHARMM --ffout CHARMM --chain tmp_pdb2pqr.pdb 00_1y0m.pqr
    Succeed to read file 00_1y0m.pqr ,  996 atoms found
    Chain: A  Residue: 0 to 60
    Succeed to save file 01_1y0m_good_his.pdb
    -Create topologie
    gmx pdb2gmx -f 01_1y0m_good_his.pdb -o 1y0m_pdb2gmx.pdb -p 1y0m_pdb2gmx.top -i \
1y0m_posre.itp -water tip3p -ff charmm36-jul2017 -ignh -vsite hydrogens
    Molecule topologie present in 1y0m_pdb2gmx.top , extract the topologie in a separate \
file: 1y0m_pdb2gmx.itp
    Protein_chain_A
    -ITP file: 1y0m_pdb2gmx.itp
    -molecules defined in the itp file:
    * Protein_chain_A
    Rewrite topologie: 1y0m_pdb2gmx.top
    >>> ######################################
    >>> ### Monitor an energy minimisation ###
    >>> ######################################
    >>> monitor = {'function': print_log_file,\
           'terms':['Potential'],\
           'file_check_ext':'log'}
    >>> prot.em(out_folder=os.path.join(TEST_OUT, 'em_SH3'), nsteps=100,\
    constraints='none', create_box_flag=True, monitor=monitor, nstlog=10)
    -Create pbc box
    gmx editconf -f .../top_SH3/1y0m_pdb2gmx.pdb -o .../top_SH3/1y0m_pdb2gmx_box.pdb -bt dodecahedron -d 1.0
    -Create the tpr file  1y0m.tpr
    gmx grompp -f 1y0m.mdp -c ../top_SH3/1y0m_pdb2gmx_box.pdb -r ../top_SH3/1y0m_pdb2gmx_box.pdb -p ../top_SH3/1y0m_pdb2gmx.top -po out_1y0m.mdp -o 1y0m.tpr -maxwarn 1
    -Launch the simulation 1y0m.tpr
    gmx mdrun -s 1y0m.tpr -deffnm 1y0m -nt 0 -ntmpi 0 -nsteps -2 -nocopyright
    Energy could not be extracted, simulation is probably finished.


    """

    log_to_check = func_input_dict['log']
    time_modif = None
    file_time = None
    count = 1
    if 'refresh_time' in func_input_dict:
        refresh_time = func_input_dict['refresh_time']
    else:
        refresh_time = 1.0


    while proc.poll() is None:

        time.sleep(refresh_time)
        count += 1

        if os_command.check_file_exist(log_to_check):
            file_time = os.stat(log_to_check).st_mtime

        if time_modif != file_time:

            time_modif = file_time

            log_dict = extract_log_dict(func_input_dict)
            try:
                print("time = {:6.1f}".format(log_dict['time']), end='')
                for keys in func_input_dict['terms']:
                    print("  {} = {:5.1f} ".format(keys, log_dict[keys]), end='')
                print()
            except KeyError:
                print('Energy could not be extracted, simulation is probably finished.')

if __name__ == "__main__":

    import doctest
    import shutil

    TEST_DIR = 'gromacs_py_test_out'
    TEST_OUT = os.path.join(TEST_DIR, 'monitor')

    def getfixture(*args):
        return TEST_OUT

    print("-Test os_command module:")

    print("monitor:  \t", doctest.testmod())
    # Erase all test files
    shutil.rmtree(TEST_DIR, ignore_errors=True)
