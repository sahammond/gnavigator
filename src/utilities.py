"""
Created on 23 May 2018

@author: S. Austin Hammond

Helper and formatting functions.

"""
import os
import subprocess
import re
import argparse

from time import localtime, strftime
from sys import argv

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


def jira_formatter(num_pct_tuple):
    # num_pct_tuple is tuple of (num, pct)
    num = num_pct_tuple[0]
    pct = num_pct_tuple[1]
    outbuff = str(num) + " " + "(" + str(pct) + "%)"
    return outbuff


def table_formatter(results_tuple):
    # results_tuple is from one of the check* functions
    # (cDNA, scaffold, status), or
    #  (cDNA, scaffold1;scaffold2..., status)
    nam = results_tuple[0]
    scaf = results_tuple[1].split(";")
    stat = results_tuple[2]

    for entry in scaf:
        outbuff = "\t".join([nam, stat, entry])
        yield outbuff


def table_formatter_wGM(results_tuple, genetic_map):
    # results_tuple is from one of the check* functions
    # (cDNA, scaffold, status),
    ##  or (cDNA, scaffold1;scaffold2..., status)
    nam = results_tuple[0]
    scaf = results_tuple[1].split(";")
    stat = results_tuple[2]
    if nam in set(genetic_map.cDNA.tolist()):
        _cdna = genetic_map[genetic_map.cDNA.isin([nam])]
        lg = str(_cdna.LG.iloc[0])
    else:
        lg = "unknownLG"

    for entry in scaf:
        outbuff = "\t".join([nam, stat, entry, lg])
        yield outbuff


def LG_table_formatter(results_tuple):
    # results_tuple is from one of the check* functions
    # (cDNA, scaffold, status),
    ##  or (cDNA, scaffold1;scaffold2..., status)
    nam = results_tuple[0]
    cDNA = " ".join(results_tuple[1].split(";"))
    stat = results_tuple[2]
    lg = results_tuple[3]

    outbuff = "\t".join([nam, cDNA, stat, lg])
    yield outbuff


def report_time():
    rep = ' '.join(["Current time:",
                    strftime("%Y-%m-%d %H:%M:%S", localtime())])
    return rep


def report_cmd():
    rep = '## Gnavigator command was: ' + ' '.join(argv)
    return rep
