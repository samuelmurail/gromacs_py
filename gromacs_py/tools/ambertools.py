#!/usr/bin/env python3

# coding: utf-8

""" Collection of function to monitor a simulation in real time.
"""

import os
import logging

from os_command_py import os_command

# Logging
logger = logging.getLogger(__name__)


# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

# Test folder path
MONITOR_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(MONITOR_LIB_DIR, "../test_files/"))

# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    logger.info("Ambertools cannot be found")
    REDUCE_BIN = ""
else:
    REDUCE_BIN = os_command.which('reduce')


def add_hydrogen(pdb_in, pdb_out, check_file_out=True, **reduce_options):
    """Add hydrogen to a pdb file using the ``reduce`` software:

    :param pdb_in: pdb input
    :type pdb_in: str

    :param pdb_out: pdb output
    :type pdb_out: str

    :param check_file_out: flag to check or not if file has already been
        created. If the file is present then the command break.
    :type check_file_out: bool, optional, default=True

    :param reduce_options: Optional arguments for ``reduce``


    :Example:

    >>> from pdb_manip_py import pdb_manip
    >>> pdb_manip.show_log()
    >>> TEST_OUT = getfixture('tmpdir')
    >>> add_hydrogen(pdb_in=TEST_PATH+'/phenol.pdb',\
    pdb_out=TEST_OUT+'/phenol_h.pdb') #doctest: +ELLIPSIS
    reduce -build -nuclear .../phenol.pdb
    >>> phenol_coor = pdb_manip.Coor(TEST_OUT+'/phenol_h.pdb')  \
#doctest: +ELLIPSIS
    Succeed to read file .../phenol_h.pdb ,  13 atoms found

    """
    # Check if output files exist:
    if check_file_out and os.path.isfile(pdb_out):
        logger.info("PDB files not created, {} already exist".format(
            pdb_out))
        return

    # Define reduce command:
    cmd_reduce = os_command.Command([REDUCE_BIN,
                                     "-build",
                                     "-nuclear",
                                     pdb_in], **reduce_options)

    cmd_reduce.display()
    return_code = cmd_reduce.run(display=False, out_data=True)

    filout = open(pdb_out, 'w')
    filout.write(return_code['stdout'])

    logger.info("Succeed to save file %s" % os.path.relpath(pdb_out))

    return


# antechamber -i phenol_h.pdb -fi pdb -o phenol_h.mol2 -fo mol2 -c bcc -s 2
