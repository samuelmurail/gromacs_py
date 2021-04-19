#!/usr/bin/env python3
# coding: utf-8

""" gmx library include the gromacs system class ``GmxSys``, as well as
topologie ``TopSys`` and itp.
"""

import logging

from .itp import Itp
from .topsys import TopSys
from .gmxsys import GmxSys, GROMACS_MOD_DIRNAME, GMX_BIN, gmx_version

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "2.0.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Prototype"

# Logging
logger = logging.getLogger(__name__)


def show_log(pdb_manip_log=True):

    from . import gmxsys
    gmxsys.show_log()
