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
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(
                sys.argv[0]))) + '/gnavigator/src')

import classify
import checklg
import reporting
import utilities as util
import config
import alignment

try:
    import pandas as pd
except ImportError:
    msg = ''.join(["""ERROR: The pandas module is required, but was not found.""",
                   """ Please install and try again."""])
    sys.exit(msg)


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
        alignment.run_gmap(args.prefix, dbDir, dbDflag, dbName, dbNflag, threads, cDNA,
                           genome, trans_mode)
        # re-check what alignments we have
        checkU, checkM, checkD = preflight(args.prefix)

    # setup results dict
    cDNA_dict = {'Complete':[], 'Duplicated':[], 'Partial':[], 'Fragmented':[],
                'Poorly mapped':[]}

    # load the fresh alignments
    uniqDat, duplDat, tlocDat = load_data(checkU, checkM, checkD, args.prefix)

    # run assessment
    print '\n=== Evaluating alignments ==='
    cDNA_res = assess(checkU, checkM, checkD, uniqDat, tlocDat, duplDat, cDNA_dict,
                      ident_thold, cov_thold)
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
    # if no genetic map data, write out the cDNA results
    else:
        # belatedly output the cDNA results without GM info
        reporting.output_cDNA(args.prefix, cDNA_res)
        reporting.report_cDNA(args.prefix, cDNA_res, TOT)


def get_args():
    """Get arguments from command line"""

    parser = argparse.ArgumentParser(
              description='Assess assembly quality & completeness using cDNA sequences')
    parser.add_argument('cDNA', help='FASTA file of cDNA sequences to align to assembly')
    parser.add_argument('genome', help='FASTA file of genome assembly to assess')
    parser.add_argument('-r', '--transcriptome', action='store_true',
                        help='Transcriptome assessment mode. See manual for details. [off]')
    parser.add_argument('-p', '--prefix', default='gnavigator',
                        help='Prefix to use for intermediate and output files [gnavigator]')
    parser.add_argument('-d', '--db_dir',
                        help='Path to directory containing prebuilt GMAP index [optional]')
    parser.add_argument('-n', '--db_name', help='Name of prebuilt GMAP index [optional]')
    parser.add_argument('-t', '--threads', action='store', default='1', type=str,
                        help='Number of threads for GMAP alignment [1]')
    parser.add_argument('-m', '--genetic_map',
                        help='Genetic map file as tsv with LG:cDNA pairs [optional]')
    parser.add_argument('-i', '--identity', action='store', type=float, default=0.95,
                        help='Minimum identity threshold [0.95]')
    parser.add_argument('-c', '--coverage', action='store', type=float, default=0.95,
                        help='Minimum coverage threshold [0.95]')

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
            print ''.join(['Gnavigator found pre-existing GMAP alignment results.',
                           ' Will use the following files:'])
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

def assess(checkU, checkD, checkM, uniqDat, tlocDat, duplDat, cDNA_results,
           ident_thold, cov_thold):
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
        if rep == 'Same LG, expected order':
            gm_res['goodLG'].append(res)
        elif rep == 'Same LG, unexpected order':
            gm_res['WO_LG'].append(res)
        elif rep == 'Different LG':
            gm_res['diffLG'].append(res)
        elif rep == 'Same LG, order undetermined':
            gm_res['undet'].append(res)

    return gm_res


if __name__ == '__main__':
    main()
