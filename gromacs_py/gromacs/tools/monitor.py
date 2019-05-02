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

def simulation_plot(proc, function_list, func_input_dict, refresh_time=1.0):
    """ This function is used for monitoring a simulation in real time.
    Function can be excecuted by the gromacs.tools.os_command.run_background() function.
    The function monitors a trajectory file, and launch the analysis if the file has been modified.
    It can plot as function of time an analysis of a simulation.
    Analysis is passed as input function.
    """

    traj_to_check = func_input_dict['xtc']
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
    #fig = plt.figure()
    fig, axarr = plt.subplots(len(function_list), sharex=True)

    for i, function in enumerate(function_list):
        y_list.append([])
        x_list.append([])
        axarr[i].set_xlabel('time (ns)')
        axarr[i].set_ylabel(function['label'])
        axarr[i].plot(x_list[i], y_list[i],
                   'ko-', markersize=2,
                   linewidth=0.5,
                   color='blue')
    # show the window
    # figure will be in foreground, but the user may move it to background
    fig.show() 
    fig.canvas.set_window_title(func_input_dict['xtc'][:-4])

    while proc.poll() is None:

        time.sleep(refresh_time)
        count += 1

        if os_command.check_file_exist(traj_to_check):
            file_time = os.stat(traj_to_check).st_mtime
        
        if time_modif != file_time:
            
            time_modif = file_time

            for i, function in enumerate(function_list):
                anal = function['name'](func_input_dict)
            
                x_list[i].append(anal['time'])
                y_list[i].append(anal['dist'])

                axarr[i].lines[0].set_data(x_list[i], y_list[i]) # set plot data
                axarr[i].relim()                  # recompute the data limits
                axarr[i].autoscale_view()         # automatic axis scaling
            fig.canvas.flush_events()   # update the plot and take care of window events (like resizing etc.)
