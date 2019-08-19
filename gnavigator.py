#!/usr/bin/env python3

"""
Created on 23 May 2018

@author: S. Austin Hammond

Use cDNA sequences +/- genetic map information to assess a genome
 assembly for completeness +/- correctness
"""

import sys
import os

# ensure gnavigator's src directory is in PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(
                sys.argv[0]))) + '/gnavigator/src')

import classify
import reporting
import utilities as util
import config
import alignment
import parseGmap
import geneticmap as gm

try:
    import pandas as pd
except ImportError:
    msg = ''.join(["""ERROR: The pandas module is required, but was not found.""",
                   """ Please install and try again."""])
    sys.exit(msg)


def main():
    args = util.get_args()

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
    checkU, checkM, checkD = util.preflight(args.prefix, 'pre')

    # load the alignments if pre-existing ones were found
    if checkU and checkM and checkD:
        uniqDat, duplDat, tlocDat = util.load_data(checkU, checkM, checkD, args.prefix)
    # if not, then proceed with alignments
    else:
        # generate the gmap config file if it doesn't exist yet
        if config.check_config():
            pass # found gmap, don't panic
        else:
            print(''.join(['ERROR: Failed to locate gmap binaries. Please install',
                           ' or specify path in gmap_config.txt manually.\n',
                            'e.g. export PATH=/home/myuser/gmap/bin:$PATH']))
            sys.exit(1)
        # run the alignments
        alignment.run_gmap(args.prefix, dbDir, dbDflag, dbName, dbNflag, threads, cDNA,
                           genome, trans_mode)
        # re-check what alignments we have
        checkU, checkM, checkD = util.preflight(args.prefix)

    # setup results dict
    cDNA_dict = {'Complete':[], 'Duplicated':[], 'Partial':[], 'Fragmented':[],
                'Poorly mapped':[]}

    # load the fresh alignments
    uniqDat, duplDat, tlocDat = util.load_data(checkU, checkM, checkD, args.prefix)

    # run assessment
    print('\n=== Evaluating alignments ===')
    cDNA_res = classify.assess(checkU, checkM, checkD, uniqDat, tlocDat, duplDat, cDNA_dict,
                      ident_thold, cov_thold)
    print('Done!')
    print(util.report_time())

    # count total number of query sequences
    check_missing = reporting.find_missing(cDNA, cDNA_res)
    TOT = check_missing[1]
    cDNA_res['Missing'] = check_missing[0]

    # load genetic map data
    if check_gm:
        if not checkU:
            print(''.join(['WARNING: There were no uniquely-aligned cDNAs',
                          ' detected, so the genetic map analysis will',
                          ' not be performed']))
            sys.exit(2)
        else:
            mapDat, uMap, uniqDatMap_select = gm.load_gm(gmfile, uniqDat, cDNA_res)
            # belatedly output the cDNA results with GM info
            reporting.output_cDNA(args.prefix, cDNA_res, mapDat)
            reporting.report_cDNA(args.prefix, cDNA_res, TOT)
        # check if there's anything to work with
        if len(uniqDatMap_select) == 0:
            print('ERROR: There are no cDNAs from the genetic map to evaluate.')
            print(''.join(['This can happen if the cDNA sequence IDs do not match those',
                  ' in the genetic map.']))
            sys.exit(2)
        else:
            gmres = gm.assess_gm(uMap, mapDat)
            reporting.output_gm(args.prefix, gmres)
            gm_cdna_stat = reporting.report_gm_cDNA(gmres, uniqDatMap_select, args.prefix) # per cDNA
            reporting.report_gm(uniqDatMap_select, gmres, gm_cdna_stat, args.prefix) # per scaffold

            # output updated genetic map
            gnavOut = f'{args.prefix}-full-cDNA-results-table.tsv'
            uniqF = f'{args.prefix}.uniq'
            duplF = f'{args.prefix}.mult'
            tlocF = f'{args.prefix}.transloc'
            parseGmap.wrapper(gnavOut, uniqF, duplF, tlocF, gmfile, args.prefix)

# if no genetic map data, write out the cDNA results
    else:
        # belatedly output the cDNA results without GM info
        reporting.output_cDNA(args.prefix, cDNA_res)
        reporting.report_cDNA(args.prefix, cDNA_res, TOT)


if __name__ == '__main__':
    main()
