#!/usr/bin/env python2

"""
Created on 23 May 2018

@author: S. Austin Hammond

Use cDNA sequences +/- genetic map information to assess a genome
 assembly for completeness +/- correctness
"""

import sys
import os
import subprocess
import re
import argparse

# ensure gnavigator's src directory is in PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/gnavigator/src')

import classify
import checklg
import reporting
import utilities as util
import config

try:
    import pandas as pd
except ImportError:
    sys.exit("""ERROR: The pandas module is required, but was not found. Please install and try again.""")


def main():
    args = get_args()

    # process the arguments    
    cDNA = args.cDNA
    genome = args.genome
    threads = args.threads
    if args.transcriptome:
        trans_mode = True
    else:
        trans_mode = False
    if args.db_dir:
        dbDflag = True
        dbDir = args.db_dir
    else:
        dbDflag = False
        dbDir = ''.join([os.getcwd(), '/', args.prefix, '-gmap-index-dir'])
    if args.db_name:
        dbNflag = True
        dbName = args.db_name
    else:
        dbNflag = False
        dbName = '-'.join([args.prefix, 'gmap-index'])
    if args.genetic_map:
        check_gm = True
        gmfile = args.genetic_map
    else:
        check_gm = False
    if args.identity < 1:
        ident_thold = args.identity
    else:
        ident_thold = float(args.identity) / 100.0
    if args.coverage < 1:
        cov_thold = args.coverage
    else:
        cov_thold = float(args.coverage) / 100.0

    # check what needs to be done
    checkU, checkM, checkD = preflight(args.prefix, 'pre')

    # load the alignments if pre-existing ones were found
    if checkU and checkM and checkD:
        uniqDat, duplDat, tlocDat = load_data(checkU, checkM, checkD, args.prefix)
    # if not, then proceed with alignments
    else:
        # generate the gmap config file if it doesn't exist yet
        if config.check_config():
            pass # found gmap, don't panic
        else:
            print ''.join(['ERROR: Failed to locate gmap binaries. Please install',
                           ' or specify path in gmap_config.txt manually.\n',
                            'e.g. export PATH=/home/myuser/gmap/bin:$PATH'])
            sys.exit(1)
        # run the alignments
        run_gmap(args.prefix, dbDir, dbDflag, dbName, dbNflag, threads, cDNA, genome, trans_mode)
        # re-check what alignments we have
        checkU, checkM, checkD = preflight(args.prefix)

    # setup results dict
    cDNA_dict = {'Complete':[], 'Duplicated':[], 'Partial':[], 'Fragmented':[],
                'Poorly mapped':[]}

    # load the fresh alignments
    uniqDat, duplDat, tlocDat = load_data(checkU, checkM, checkD, args.prefix)

    # run assessment
    print '\n=== Evaluating alignments ==='
    cDNA_res = assess(checkU, checkM, checkD, uniqDat, tlocDat, duplDat, cDNA_dict, ident_thold, cov_thold)
    print 'Done!'
    print util.report_time()

    # count total number of query sequences
    check_missing = reporting.find_missing(cDNA, cDNA_res)
    TOT = check_missing[1]
    cDNA_res['Missing'] = check_missing[0]

    # load genetic map data
    if check_gm:
        if not checkU:
            print ''.join(['WARNING: There were no uniquely-aligned cDNAs',
                          ' detected, so the genetic map analysis will',
                          ' not be performed'])
            sys.exit(2)
        else:
            mapDat, uMap, uniqDatMap = load_gm(gmfile, uniqDat)
            # belatedly output the cDNA results with GM info
            reporting.output_cDNA(args.prefix, cDNA_res, mapDat)
            reporting.report_cDNA(args.prefix, cDNA_res, TOT)
        # check if there's anything to work with
        if len(uniqDatMap) == 0:
            print 'ERROR: There are no cDNAs from the genetic map to evaluate.'
            print ''.join(['This can happen if the cDNA sequence IDs do not match those',
                  ' in the genetic map.'])
            sys.exit(2)
        else:
            gmres = assess_gm(uMap, mapDat)
            reporting.output_gm(args.prefix, gmres)
            reporting.report_gm(uMap, gmres, args.prefix)
    else:
        # belatedly output the cDNA results without GM info
        reporting.output_cDNA(args.prefix, cDNA_res)
        reporting.report_cDNA(args.prefix, cDNA_res, TOT)


def get_args():
    """Get arguments from command line"""

    parser = argparse.ArgumentParser(description='Assess assembly quality and completeness using cDNA sequences')
    parser.add_argument('cDNA', help='FASTA file of cDNA sequences to align to assembly') # cDNA sequence fasta
    parser.add_argument('genome', help='FASTA file of genome assembly to assess')
    parser.add_argument('-r', '--transcriptome', help='Transcriptome assessment mode. See manual for details. [off]', action='store_true') # use nonspliced alignments
    parser.add_argument('-p', '--prefix', help='Prefix to use for intermediate and output files [gnavigator]', default='gnavigator') # prefix
    parser.add_argument('-d', '--db_dir', help='Path to directory containing prebuilt GMAP index [optional]') # gmap db dir
    parser.add_argument('-n', '--db_name', help='Name of prebuilt GMAP index [optional]') # gmap db name
    parser.add_argument('-t', '--threads', help='Number of threads for GMAP alignment [1]', action='store', default='1', type=str)
    parser.add_argument('-m', '--genetic_map', help='Genetic map file as tsv with LG:cDNA pairs [optional]')
    parser.add_argument('-i', '--identity', help='Minimum identity threshold [0.95]', action='store', default=0.95, type=float)
    parser.add_argument('-c', '--coverage', help='Minimum coverage threshold [0.95]', action='store', default=0.95, type=float)

    arguments = parser.parse_args()

    return arguments


def preflight(prefix, stage='post'):
    """check which alignment files were produced"""
    checkU = os.path.isfile(''.join([os.getcwd(), '/', prefix, '.uniq']))
    checkM = os.path.isfile(''.join([os.getcwd(), '/', prefix, '.mult']))
    checkD = os.path.isfile(''.join([os.getcwd(), '/', prefix, '.transloc']))

    if stage == 'pre':
        if checkU or checkM or checkD:
            print '\n=== Skipping GMAP alignment stage ==='
            print 'Gnavigator found pre-existing GMAP alignment results. Will use the following files:'
            if checkU:
                print ''.join([os.getcwd(), '/', prefix, '.uniq'])
            if checkM:
                print ''.join([os.getcwd(), '/', prefix, '.mult'])
            if checkD:
                print ''.join([os.getcwd(), '/', prefix, '.transloc'])
            print util.report_time()
    elif stage == 'post':
        pass

    return (checkU, checkM, checkD)


def run_gmap(prefix, dbDir, dbDflag, dbName, dbNflag, threads, cDNA, genome, trans_mode):
    gnavigator_path = re.sub('gnavigator.py', '', os.path.realpath(__file__))

    # detect pre-existing index, use this check later (ref153positions is final index file created)
    checkI = os.path.isfile(''.join([os.getcwd(), '/', prefix, '-gmap-index-dir/',
                                     prefix, '-gmap-index/', prefix, '-gmap-index.ref153positions']))
    # make log file names
    indexlog = '-'.join([prefix, 'gmap', 'index.log'])
    alignlog = '-'.join([prefix, 'gmap', 'alignment.log'])
    # check if user supplied an index
    if dbDflag and dbNflag:
        print '\n=== Skipping GMAP index construction ==='
        print 'Gnavigator will use the user-specified index:'
        print dbDir
        print util.report_time()
    # if not, check if index made already
    elif checkI:
        print '\n=== Skipping GMAP index construction ==='
        print 'Gnavigator found a pre-existing GMAP index:'
        print ''.join([os.getcwd(), '/', prefix, '-gmap-index-dir'])
        print util.report_time()      
    # otherwise, make gmap index
    else:
        print '\n=== Building GMAP database ==='
        print util.report_time()
        try:
            index_cmd = [gnavigator_path + 'bin/build-index.sh', dbDir, dbName, genome, indexlog]
            subprocess.check_call(index_cmd)
        except subprocess.CalledProcessError:
            print '\nERROR: Failed to build GMAP index.'
            print 'Make sure that the genome file exists.'
            sys.exit(1)
        print 'Done!'
    # run gmap alignment
    print '\n=== Performing GMAP alignments ==='
    print util.report_time()
    if trans_mode:    
        print 'Running in transcriptome assessment mode. Will run GMAP without splicing.'
        try:
            # try running regular gmap first; will fail if genome too big
            gmap_cmd = [gnavigator_path + 'bin/run-gmap.sh', dbDir, dbName, threads, prefix, cDNA, alignlog, 'N']
            subprocess.check_call(gmap_cmd)
        except subprocess.CalledProcessError:
            try:
                # genome was probably too big. Use gmapl
                gmapl_cmd = [gnavigator_path + 'bin/run-gmapl.sh', dbDir, dbName, threads, prefix, cDNA, alignlog, 'N']
                subprocess.check_call(gmapl_cmd)
            except subprocess.CalledProcessError as e:
                print '\nERROR: Failed to perform GMAP alignment.'
                print 'Make sure that the cDNA file exists.'
                print 'GMAP error text:'
                print e.output
                sys.exit(1)
    else:
        # not in transcriptome mode; run WITH splicing
        try:
            # try running regular gmap first; will fail if genome too big
            gmap_cmd = [gnavigator_path + 'bin/run-gmap.sh', dbDir, dbName, threads, prefix, cDNA, alignlog]
            subprocess.check_call(gmap_cmd)
        except subprocess.CalledProcessError:
            try:
                # genome was probably too big. Use gmapl
                gmapl_cmd = [gnavigator_path + 'bin/run-gmapl.sh', dbDir, dbName, threads, prefix, cDNA, alignlog]
                subprocess.check_call(gmapl_cmd)
            except subprocess.CalledProcessError as e:
                print '\nERROR: Failed to perform GMAP alignment.'
                print 'Make sure that the cDNA file exists.'
                print 'GMAP error text:'
                print e.output
                sys.exit(1)
    print 'Done!'


def load_data(checkU, checkM, checkD, prefix):
    """read in the data and define extent"""
    col_names = ['matches', 'mismatches', 'repmatches', 'ncount',
                 'qnuminsert', 'qbaseinsert', 'tnuminsert', 'tbaseinsert',
                 'strand', 'qname', 'qsize', 'qstart', 'qend', 'tname',
                 'tsize', 'tstart', 'tend', 'blockcount', 'blocksizes',
                 'qstarts', 'tstarts']
    if checkU:
        uniqDat = pd.read_csv('.'.join([prefix, 'uniq']), sep='\t',
                              comment='#', low_memory=False, header=None,
                              names=col_names)
    if checkM:
        duplDat = pd.read_csv('.'.join([prefix, 'mult']), sep='\t',
                              comment='#', low_memory=False, header=None,
                              names=col_names)
    if checkD:
        tlocDat = pd.read_csv('.'.join([prefix, 'transloc']), sep='\t',
                              comment='#', low_memory=False, header=None,
                              names=col_names)

    return (uniqDat, duplDat, tlocDat)

def assess(checkU, checkD, checkM, uniqDat, tlocDat, duplDat, cDNA_results, ident_thold, cov_thold):
    # apply check_complete to whole set
    if checkU:
        for rec in uniqDat.itertuples():
            res = classify.check_aln(rec, 'assess', ident_thold, cov_thold)
            cDNA_results[res[2]].append(res) # append results tuple

    # apply check_frag to whole set
    if checkD:
        for qry in tlocDat.qname.unique():
            this_qry = tlocDat[tlocDat.qname == qry]
            frags = []
            for rec in this_qry.itertuples():
                frags.append(rec)
            res = classify.check_frag(frags, ident_thold, cov_thold)
            cDNA_results[res[2]].append(res)

    # apply check_dupl to whole set
    if checkM:
        for qry in duplDat.qname.unique():
            this_qry = duplDat[duplDat.qname == qry]
            frags = []
            for rec in this_qry.itertuples():
                frags.append(rec)
            res = classify.check_dupl(frags, ident_thold, cov_thold)
            cDNA_results[res[2]].append(res)

    return cDNA_results


def load_gm(gmfile, uniqDat):
    """Read in genetic map. Expect format LG\tcM\tcDNA"""
    mapDat = pd.read_csv(gmfile, sep='\t', comment='#', low_memory=False,
                         header=None, names=['LG', 'cM', 'cDNA'])
    mapDat = mapDat.sort_values(['LG', 'cM'], ascending=True)
    # limit genetic map analysis to complete (i.e. single) cDNAs
    map_cDNA = set(mapDat.cDNA.tolist())
    uniqDatMap = uniqDat[uniqDat.qname.isin(map_cDNA)]
    uMap = uniqDatMap[uniqDatMap.tname.duplicated(keep=False)]

    return (mapDat, uMap, uniqDatMap)


def assess_gm(uMap, mapDat):
    """apply check_LG to whole uniq set"""
    gm_res = {'goodLG':[], 'WO_LG':[], 'diffLG':[], 'undet':[]}
    for rec in uMap.tname.unique():
        thisScaf = uMap[uMap.tname.isin([rec])]
        res = checklg.check_LG(thisScaf, mapDat)
        rep = res[2]
        if rep == 'Same LG, right order':
            gm_res['goodLG'].append(res)
        elif rep == 'Same LG, wrong order':
            gm_res['WO_LG'].append(res)
        elif rep == 'Different LG':
            gm_res['diffLG'].append(res)
        elif rep == 'Same LG, order undetermined':
            gm_res['undet'].append(res)

    return gm_res


if __name__ == '__main__':
    main()
