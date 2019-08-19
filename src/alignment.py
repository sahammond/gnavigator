"""
Created on 6 June 2018

@author: S. Austin Hammond

Functions to run the alignments.
"""

import sys
import re
import os
import subprocess

# ensure gnavigator's src directory is in PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(
                sys.argv[0]))) + '/gnavigator/src')

import utilities as util

def run_gmap(prefix, dbDir, dbDflag, dbName, dbNflag, threads, cDNA, genome, trans_mode):
    gnavigator_path = re.sub('src/alignment.pyc', '', os.path.realpath(__file__))
    gnavigator_path = re.sub('src/alignment.py', '', gnavigator_path)
    # detect pre-existing index, use this check later
    # n.b. ref153positions is final index file created
    checkI = os.path.isfile(''.join([os.getcwd(), '/', prefix, '-gmap-index-dir/',
                                     prefix, '-gmap-index/', prefix, 
                                     '-gmap-index.ref153positions']))
    # make log file names
    indexlog = '-'.join([prefix, 'gmap', 'index.log'])
    alignlog = '-'.join([prefix, 'gmap', 'alignment.log'])
    # check if user supplied an index
    if dbDflag and dbNflag:
        print('\n=== Skipping GMAP index construction ===')
        print('Gnavigator will use the user-specified index:')
        print(dbDir)
        print(util.report_time())
    # if not, check if index made already
    elif checkI:
        print('\n=== Skipping GMAP index construction ===')
        print('Gnavigator found a pre-existing GMAP index:')
        print(''.join([os.getcwd(), '/', prefix, '-gmap-index-dir']))
        print(util.report_time())
    # otherwise, make gmap index
    else:
        print('\n=== Building GMAP database ===')
        print(util.report_time())
        try:
            index_cmd = [gnavigator_path + 'bin/build-index.sh', dbDir,
                         dbName, genome, indexlog]
            subprocess.check_call(index_cmd)
        except subprocess.CalledProcessError:
            print('\nERROR: Failed to build GMAP index.')
            print('Make sure that the genome file exists.')
            sys.exit(1)
        print('Done!')
    # run gmap alignment
    print('\n=== Performing GMAP alignments ===')
    print(util.report_time())
    if trans_mode:    
        print('Running in transcriptome assessment mode. Will run GMAP without splicing.')
        try:
            # try running regular gmap first; will fail if genome too big
            gmap_cmd = [gnavigator_path + 'bin/run-gmap.sh', dbDir, dbName,
                        threads, prefix, cDNA, alignlog, 'N']
            subprocess.check_call(gmap_cmd)
        except subprocess.CalledProcessError:
            try:
                # genome was probably too big. Use gmapl
                gmapl_cmd = [gnavigator_path + 'bin/run-gmapl.sh', dbDir, dbName,
                             threads, prefix, cDNA, alignlog, 'N']
                subprocess.check_call(gmapl_cmd)
            except subprocess.CalledProcessError as e:
                print('\nERROR: Failed to perform GMAP alignment.')
                print('Make sure that the cDNA file exists.')
                print('GMAP error text:')
                print(e.output)
                sys.exit(1)
    else:
        # not in transcriptome mode; run WITH splicing
        try:
            # try running regular gmap first; will fail if genome too big
            gmap_cmd = [gnavigator_path + 'bin/run-gmap.sh', dbDir, dbName, threads,
                        prefix, cDNA, alignlog]
            subprocess.check_call(gmap_cmd)
        except subprocess.CalledProcessError:
            try:
                # genome was probably too big. Use gmapl
                gmapl_cmd = [gnavigator_path + 'bin/run-gmapl.sh', dbDir, dbName,
                             threads, prefix, cDNA, alignlog]
                subprocess.check_call(gmapl_cmd)
            except subprocess.CalledProcessError as e:
                print('\nERROR: Failed to perform GMAP alignment.')
                print('Make sure that the cDNA file exists.')
                print('GMAP error text:')
                print(e.output)
                sys.exit(1)
    print('Done!')
