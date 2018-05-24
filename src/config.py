#!/usr/bin/env python2

"""
Created on 23 May 2018

@author: S. Austin Hammond

Generate an initial configuration script for GMAP
"""

import os
import sys


def check_config():
    """create a config file if none exists"""
    conf_file = ''.join([os.getcwd(), '/gmap_config.txt'])
    if not os.path.isfile(conf_file):
        make_config()


def make_config():
    """find default gmap binaries"""
    gmap = sys.call(['which', 'gmap'])
    gmap_dir = gmap[:-5]

    outname = 'gmap_config.txt'

    with open(outname, 'w') as outfile:
        path_str = ''.join(['export ', 'PATH=', gmap_dir, ':$PATH'])
        print >> outfile, path_str
