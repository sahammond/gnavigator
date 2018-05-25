#!/usr/bin/env python2

"""
Created on 23 May 2018

@author: S. Austin Hammond

Generate an initial configuration script for GMAP
"""

import os
import subprocess
import re


def check_config():
    """create a config file if none exists"""
    gnavigator_path = re.sub('src/config.pyc', '', os.path.realpath(__file__))
    gnavigator_path = re.sub('src/config.py', '', gnavigator_path)
    conf_file = ''.join([gnavigator_path, 'gmap_config.txt'])
    if not os.path.isfile(conf_file):
        msg = ' '.join([conf_file, ' not detected. Auto-generating'])
        print msg
        status = make_config()
        if status:
            msg = ' '.join(['Created', conf_file])
            print msg
    else:
        msg = ' '.join(['Found', conf_file])
        print msg
        status = True

    return status


def make_config():
    """find default gmap binaries"""
    gmap = subprocess.Popen(['which', 'gmap'], stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
    out, err = gmap.communicate()

    if not err:
        pass_flag = True
        gmap_dir = re.sub('/gmap$', '', out).strip('\n')
    else:    
        pass_flag = False

    if pass_flag:
        gnavigator_path = re.sub('src/config.pyc', '', os.path.realpath(__file__))
        gnavigator_path = re.sub('src/config.py', '', gnavigator_path)
        outname = gnavigator_path + 'gmap_config.txt'
        with open(outname, 'w') as outfile:
            path_str = ''.join(['export ', 'PATH=', gmap_dir, ':$PATH'])
            print >> outfile, path_str

    # check flag in main script and die if necessary
    return pass_flag
