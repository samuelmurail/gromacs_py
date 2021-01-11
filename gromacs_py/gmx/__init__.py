#!/usr/bin/env python3
# coding: utf-8

""" gmx library include the gromacs system class ``GmxSys``, as well as
topologie ``TopSys`` and itp.
"""

import sys
import logging

from .itp import Itp
from .topsys import TopSys
from .gmxsys import GmxSys, GROMACS_MOD_DIRNAME

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "1.2.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Prototype"

# Logging
logger = logging.getLogger(__name__)


def show_log(pdb_manip_log=True):
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.INFO)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))

    # Show log of top sys:
    from . import gmxsys
    gmxsys.show_log()
