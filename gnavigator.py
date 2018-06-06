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
            output_cDNA(args.prefix, cDNA_res, mapDat)
            report_cDNA(args.prefix, cDNA_res, TOT)
        # check if there's anything to work with
        if len(uniqDatMap) == 0:
            print 'ERROR: There are no cDNAs from the genetic map to evaluate.'
            print ''.join(['This can happen if the cDNA sequence IDs do not match those',
                  ' in the genetic map.'])
            sys.exit(2)
        else:
            gmres = assess_gm(uMap, mapDat)
            output_gm(args.prefix, gmres)
            report_gm(uMap, gmres, args.prefix)
    else:
        # belatedly output the cDNA results without GM info
        output_cDNA(args.prefix, cDNA_res)
        report_cDNA(args.prefix, cDNA_res, TOT)


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


def output_cDNA(prefix, cDNA_results, gm=''):
    # write out cDNA:scaffold mappings
    full_out = '-'.join([prefix, 'full-cDNA-results-table.tsv'])
    with open(full_out, 'w') as outfile:
        if len(gm) > 0:
            header = '\t'.join(['# cDNA ID', 'Status', 'Scaffold',
                                'Linkage group'])
            print >> outfile, header
            for status, result in cDNA_results.items():
                for res in result:
                    for t in util.table_formatter_wGM(res, gm):
                        print >> outfile, t
        else:
            header = '\t'.join(['# cDNA ID', 'Status', 'Scaffold'])
            print >> outfile, header
            for status, result in cDNA_results.items():
                for res in result:
                    for t in util.table_formatter(res):
                        print >> outfile, t


def report_cDNA(prefix, cDNA_results, TOT):
    """calc percentages and report cDNA results"""
    num_complete = len(cDNA_results['Complete'])
    num_duplicated = len(cDNA_results['Duplicated'])
    num_partial = len(cDNA_results['Partial'])
    num_fragmented = len(cDNA_results['Fragmented'])
    num_poor = len(cDNA_results['Poorly mapped'])
    num_missing = len(cDNA_results['Missing'])
    rate_complete = float(num_complete) / float(TOT)
    rate_duplicated = float(num_duplicated) / float(TOT)
    rate_partial = float(num_partial) / float(TOT)
    rate_fragmented = float(num_fragmented) / float(TOT)
    rate_poor = float(num_poor) / float(TOT)
    rate_missing = float(num_missing) / float(TOT)
    pct_complete = round(100.0 * rate_complete, 2)
    pct_duplicated = round(100.0 * rate_duplicated, 2)
    pct_partial = round(100.0 * rate_partial, 2)
    pct_fragmented = round(100 * rate_fragmented, 2)
    pct_poor = round(100 * rate_poor, 2)
    pct_missing = round(100 * rate_missing, 2)

    # report if the right number of sequences have a result
    num_counted = sum([num_complete, num_duplicated, num_fragmented,
                       num_partial, num_poor, num_missing])
    rate_counted = float(num_counted) / float(TOT)
    pct_counted = round(100 * rate_counted, 2)

    # write to tsv
    tsvout = '-'.join([prefix, 'results.tsv'])
    with open(tsvout, 'w') as outfile:
        header = '\t'.join(['', 'Complete', 'Duplicated', 'Fragmented',
                            'Partial', 'Poorly Mapped', 'Missing',
                            'Total cDNAs searched'])
        nums = '\t'.join([str(x) for x in ['Number', num_complete, num_duplicated,
                                           num_fragmented, num_partial, num_poor,
                                           num_missing, num_counted]])
        pcts = '\t'.join([str(x) for x in ['Percent', pct_complete, pct_duplicated,
                                           pct_fragmented, pct_partial, pct_poor,
                                           pct_missing, pct_counted]])
        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts
    jiraout = '-'.join([prefix, 'results.jira'])
    with open(jiraout, 'w') as outfile:
        header = '||'.join(['', 'Complete', 'Duplicated', 'Fragmented',
                            'Partial', 'Poorly Mapped', 'Missing',
                            'Total cDNAs searched', ''])
        nums = [num_complete, num_duplicated, num_fragmented, num_partial, num_poor,
                num_missing, num_counted]
        pcts = [pct_complete, pct_duplicated, pct_fragmented, pct_partial, pct_poor,
                pct_missing, pct_counted]
        res = '|' + '|'.join([util.jira_formatter(x) for x in zip(nums, pcts)]) + '|'
        
        print >> outfile, header
        print >> outfile, res

    # print to STDOUT
    print '\n=== GNAVIGATOR cDNA RESULTS ==='
    print '%s (%s%%) complete sequences' % (num_complete, pct_complete)
    print '%s (%s%%) duplicated sequences' % (num_duplicated, pct_duplicated)
    print '%s (%s%%) fragmented sequences' % (num_fragmented, pct_fragmented)
    print '%s (%s%%) partial sequences' % (num_partial, pct_partial)
    print '%s (%s%%) poorly mapped sequences' % (num_poor, pct_poor)
    print '%s (%s%%) missing sequences' % (num_missing, pct_missing)
    print '%s (%s%%) sequences were evaluated' % (num_counted, pct_counted)
    print util.report_time()


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


# write out cDNA:scaffold mappings
def output_gm(prefix, gm_res):
    header = '\t'.join(['# Scaffold', 'cDNA IDs', 'Status', 'Linkage group(s)'])
    full_out = '-'.join([prefix, 'full-genetic-map-results-table.tsv'])
    with open(full_out, 'w') as outfile:
        print >> outfile, header
        for status, result in gm_res.items():
            for res in result:
                for t in util.LG_table_formatter(res):
                    print >> outfile, t


def report_gm(uMap, gm_results, prefix):
    """report summary of genetic map results"""
    num_scaff_toCheck = len(uMap.tname.unique())
    num_goodLG = len(gm_results['goodLG'])
    num_WO_LG = len(gm_results['WO_LG'])
    num_diffLG = len(gm_results['diffLG'])
    num_undet = len(gm_results['undet'])
    num_scaff_checked = num_goodLG + num_WO_LG + num_diffLG + num_undet
    if num_scaff_toCheck == num_scaff_checked:
        rate_LGscaff = float(num_scaff_checked) / float(num_scaff_toCheck)
        rate_goodLG = float(num_goodLG) / float(num_scaff_checked)
        rate_WO_LG = float(num_WO_LG) / float(num_scaff_checked)
        rate_diffLG = float(num_diffLG) / float(num_scaff_checked)
        rate_undet = float(num_undet) / float(num_scaff_checked)
        pct_LGscaff = round(100.0 * rate_LGscaff, 2) 
        pct_goodLG = round(100.0 * rate_goodLG, 2)
        pct_WO_LG = round(100.0 * rate_WO_LG, 2)
        pct_diffLG = round(100.0 * rate_diffLG, 2)
        pct_undet = round(100.0 * rate_undet, 2)
    else:
        print ''.join('Not all scaffolds to be checked against genetic map'
                      ' were successfully checked.')
        print 'Maybe something is wrong with the input data?'
        print util.report_time()
        sys.exit(2)
    
    # write to tsv
    tsvout = '-'.join([prefix, 'genetic-map-results.tsv'])
    with open(tsvout, 'w') as outfile:
        header = '\t'.join(['', 'Same LG, right order', 'Same LG, wrong order',
                            'Different LG', 'Same LG, undetermined order',
                            'Total scaffolds analyzed'])
        nums = '\t'.join([str(x) for x in ['Number', num_goodLG, num_WO_LG,
                                           num_diffLG, num_undet, num_scaff_checked]])
        pcts = '\t'.join([str(x) for x in ['Percent', pct_goodLG, pct_WO_LG,
                                           pct_diffLG, pct_undet, pct_LGscaff]])
        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts
    jiraout = '-'.join([prefix, 'genetic-map-results.jira'])
    with open(jiraout, 'w') as outfile:
        header = '||'.join(['', 'Same LG, right order', 'Same LG, wrong order',
                            'Different LG', 'Same LG, undetermined order',
                            'Total scaffolds analyzed', ''])
        nums = [num_goodLG, num_WO_LG, num_diffLG, num_undet, num_scaff_checked]
        pcts = [pct_goodLG, pct_WO_LG, pct_diffLG, pct_undet, pct_LGscaff]
        res = '|' + '|'.join([util.jira_formatter(x) for x in zip(nums, pcts)]) + '|'
        print >> outfile, header
        print >> outfile, res
    print '\n=== GNAVIGATOR GENETIC MAP RESULTS ==='
    print '%s (%s%%) scaffolds had 2+ complete cDNAs from the genetic map aligned to them.' % (num_scaff_checked, pct_LGscaff)
    print '%s (%s%%) case(s) were from the same linkage group and in the expected order.' % (num_goodLG, pct_goodLG)
    print '%s (%s%%) case(s) were from the same linkage group, but NOT in the expected order.' % (num_WO_LG, pct_WO_LG)
    print '%s (%s%%) case(s) were from different linkage groups.' % (num_diffLG, pct_diffLG)
    print '%s (%s%%) case(s) were from the same linkage group but their order could not be determined.' % (num_undet, pct_undet)
    print util.report_time()


if __name__ == '__main__':
    main()
