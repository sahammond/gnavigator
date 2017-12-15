#!/usr/bin/env python

# purpose: create a tool similar to CEGMA/BUSCO that uses genetic map information to assess 
#    completeness and quality of an assembly and reporting a simple set of metrics

# want to distinguish complete, partial, fragmented, duplicated, poorly mapped, and missing
# complete: 90% of sequence is aligned with 95% identity
# partial: less than 90% of sequence is aligned to a single scaffold, other 10% unaligned
# fragmented: 90% of a sequence is aligned, but over multiple scaffolds
## would need to use check_aln in 'report' mode and assess separately
## alignment would be reported in "*.transloc"
# duplicated: as complete, but at multiple locations
## alignment would be reported in "*.mult"
# poorly mapped: as complete, partial, or fragmented, but with less than 90% identity
# missing: not aligned by gmap
## assessed separately

import sys
import os
import subprocess
import re
import argparse

import pandas as pd

# TODOs
# report which cDNAs fall into which category
# return the sequence from the assembly that corresponds to the cDNA
## related: return cDNA:scaffold relationships
# report if any genetic map cDNAs should be found between cDNAs that are observed

def check_aln(aln, mode):
    """check the alignment of a gcat map sequence"""

    matches = float(aln.matches)
    mismatches = float(aln.mismatches)
    qinserts = float(aln.qbaseinsert)
    qsize = float(aln.qsize)
    qstart = float(aln.qstart)
    qend = float(aln.qend)

    
    seg = qend - qstart
    pid = matches / seg
    # want to penalize insertions b/c reflects correctness of assembly
    pcov = (matches + mismatches - qinserts) / qsize
    
    if mode == 'assess':
        if pid >= 0.95 and pcov >= 0.95:
            return 'complete'
        elif pid >= 0.95 and pcov < 0.95:
            if pcov >= 0.5:
                return 'partial'
            else:
                return 'poorly mapped'
        else:
            return 'poorly mapped'
    
    elif mode == 'report':
        goodb = matches
        # want to penalize insertions b/c reflects correctness of assembly
        covb = matches + mismatches - qinserts
        
        return (goodb, covb, seg, qsize)


def check_frag(alns):
    """check if 'chimeric' alignment is good or not"""
    # take in a series of blat alignments and evaluate them together
    # TODO extend to assess non-chimeric alignments
    ## some cDNAs may be fragmented along >2 scaffolds, with some missing bits in the middle
    goodb = 0
    covb = 0
    seg = 0
    qsize = 0
    
    for aln in alns:
        this_aln = check_aln(aln, 'report')
        goodb += this_aln[0]
        covb += this_aln[1]
        seg += this_aln[2]
        qsize = this_aln[3]
    
    pid = goodb / seg
    pcov = covb / qsize
    
    if pid >= 0.95 and pcov >= 0.95:
        return 'fragmented'
    elif pid >= 0.95 and pcov < 0.95:
        if pcov > 0.5:
            return 'partial' # distinguish from partials in a single piece?
        else:
            return 'poorly mapped'
    else:
        return 'poorly mapped'


def check_dupl(alns):
    """check if a gcat sequence should be considered duplicated or not based on its multiple alignments"""
    # more than one alignment must be considered complete in order to be duplicated
    # one or more partial, fragmented, poorly mapped will not count
    results = []
    for aln in alns:
        this_aln = check_aln(aln, 'assess')
        results.append(this_aln)
    
    if 'complete' in results:
        num_complete = len([x for x in results if x == 'complete'])
        if num_complete > 1:
            return 'duplicated'
        elif num_complete == 1:
            return 'complete'
    else:
        if 'partial' in results:
            return 'partial'
        else:
            return 'poorly mapped'


def count_aligned(uniqA, duplA, tlocA):
    """count how many seqs have any alignment whatsoever"""
    aligned = set()
    for rec in uniqA.qname.unique():
        aligned.add(rec)
    for rec in duplA.qname.unique():
        aligned.add(rec)
    for rec in tlocA.qname.unique():
        aligned.add(rec)
        
    return len(aligned)


def check_LG(query, genetic_map):
    """check if cDNAs are from the same LG"""
    # query is a pandas array of a single scaffold's alignments
    # genetic_map is a pandas array of LG\tcM\tcDNA
    # assume for now that cDNA is in genetic_map$cDNA
    # returns 'same LG, right order', 'same LG, wrong order', or 'different LG'
    refs = query.qname.tolist()
    thisMap = genetic_map[genetic_map.cDNA.isin(refs)]
    numLG = len(thisMap.LG.unique())
    if len(thisMap) == 2:
        if numLG == 1:
            return 'same LG, right order'
        else:
            return 'different LG'
    else:
        if numLG == 1:
            # comparing two lists of the same length will return True if order is the same
            # compare both forward and reverse orders
            fwdA = query.sort_values(['tstart'])
            fwdL = fwdA.qname.tolist()
            revA = query.sort_values(['tstart'], ascending=False)
            revL = revA.qname.tolist()
            mapL = thisMap.cDNA.tolist()
            if mapL == fwdL:
                return 'same LG, right order'
            elif mapL == revL:
                return 'same LG, right order'
            else:
                return 'same LG, wrong order'
        else:
            return 'different LG'


def jira_formatter(num_pct_tuple):
    # num_pct_tuple is tuple of (num, pct)
    num = num_pct_tuple[0]
    pct = num_pct_tuple[1]
    outbuff = str(num) + " " + "(" + str(pct) + "%)"
    return outbuff


# setup parser
parser = argparse.ArgumentParser(description='Assess assembly quality and completeness using cDNA sequences')
parser.add_argument('cDNA', help='FASTA file of cDNA sequences to align to assembly') # cDNA sequence fasta
parser.add_argument('genome', help='FASTA file of genome assembly to assess')
parser.add_argument('-p', '--prefix', help='Prefix to use for intermediate and output files [gnavigator]', default='gnavigator') # prefix
parser.add_argument('-d', '--db_dir', help='Path to directory containing prebuilt GMAP index [optional]') # gmap db dir
parser.add_argument('-n', '--db_name', help='Name of prebuilt GMAP index [optional]') # gmap db name
parser.add_argument('-t', '--threads', help='Number of threads for GMAP alignment [1]', action='store', default='1')
parser.add_argument('-m', '--genetic_map', help='Genetic map file as tsv with LG:cDNA pairs [optional]')


# process args
args = parser.parse_args()
cDNA = args.cDNA
genome = args.genome
prefix = args.prefix
threads = args.threads

if args.db_dir:
    dbDir = args.db_dir
else:
    dbDir = ''.join([os.getcwd(), '/', prefix, '-gmap-index-dir'])
if args.db_name:
    dbName = args.db_name
else:
    dbName = '-'.join([prefix, 'gmap-index'])
if args.genetic_map:
    genetMap = args.genetic_map

# get path to this script, and assume that the gmap sh scripts are there too
gnavigator_path = re.sub('gnavigator.py', '', os.path.realpath(__file__))

# check if alignments have already been done
checkU = os.path.isfile(''.join([os.getcwd(), '/', prefix, ".uniq"]))
checkM = os.path.isfile(''.join([os.getcwd(), '/', prefix, ".mult"]))
checkD = os.path.isfile(''.join([os.getcwd(), '/', prefix, ".transloc"]))
if checkU and checkM and checkD:
    print "\n=== Skipping GMAP alignment stage ==="
    print "Gnavigator found pre-existing GMAP alignment results. Will use the following files:"
    print ''.join([os.getcwd(), '/', prefix, ".uniq"])
    print ''.join([os.getcwd(), '/', prefix, ".mult"])
    print ''.join([os.getcwd(), '/', prefix, ".transloc"])

# make gmap index if not supplied by user
else:
    if not args.db_dir and not args.db_name:
        print "\n ===Building GMAP database=== "
        try:
            subprocess.call([gnavigator_path + '/build-index.sh', dbDir, dbName, genome])
        except:
            print 'Failed to build GMAP index.'
            print 'Make sure that build-index.sh is in the same directory as gnavigator.'
            sys.exit(1)

# run gmap alignment
    print "\n=== Performing GMAP alignments ==="
    try:
        subprocess.call([gnavigator_path + '/run-gmap.sh', dbDir, dbName, threads, prefix, cDNA])
    except:
        print 'Failed to perform GMAP alignment.'
        print 'Make sure that run-gmap.sh is in the same directory as gnavigator.'
        sys.exit(1)

# read in the data and define extent
col_names = ['matches', 'mismatches', 'repmatches', 'ncount', 'qnuminsert', 'qbaseinsert', 'tnuminsert', 
             'tbaseinsert', 'strand', 'qname', 'qsize', 'qstart', 'qend', 'tname', 'tsize', 'tstart', 'tend',
             'blockcount', 'blocksizes', 'qstarts', 'tstarts']

uniqDat = pd.read_csv('.'.join([prefix, 'uniq']), sep='\t', comment='#', low_memory=False, header=None, names=col_names)
duplDat = pd.read_csv('.'.join([prefix, 'mult']), sep='\t', comment='#', low_memory=False, header=None, names=col_names)
tlocDat = pd.read_csv('.'.join([prefix, 'transloc']), sep='\t', comment='#', low_memory=False, header=None, names=col_names)

uniqDat['qname'] = uniqDat['qname'].str[:-2]
duplDat['qname'] = duplDat['qname'].str[:-2]
tlocDat['qname'] = tlocDat['qname'].str[:-2]

# total number of query sequences 
TOT = 27143 #TODO derive from input

### FOR DEV ###
# read in the data and define extent
#col_names = ['matches', 'mismatches', 'repmatches', 'ncount', 'qnuminsert', 'qbaseinsert', 'tnuminsert', 
#             'tbaseinsert', 'strand', 'qname', 'qsize', 'qstart', 'qend', 'tname', 'tsize', 'tstart', 'tend',
#             'blockcount', 'blocksizes', 'qstarts', 'tstarts']

#uniqDat = pd.read_csv('gmap-tests/sitka-new-strategy-postLINKS-95.uniq', sep='\t', comment='#', low_memory=False, header=None, names=col_names)
#duplDat = pd.read_csv('gmap-tests/sitka-new-strategy-postLINKS-95.mult', sep='\t', comment='#', low_memory=False, header=None, names=col_names)
#tlocDat = pd.read_csv('gmap-tests/sitka-new-strategy-postLINKS-95.transloc', sep='\t', comment='#', low_memory=False, header=None, names=col_names)

# remove the ".1" etc. added by GenBank to the query names
#uniqDat['qname'] = uniqDat['qname'].str[:-2]
#duplDat['qname'] = duplDat['qname'].str[:-2]
#tlocDat['qname'] = tlocDat['qname'].str[:-2]

# total number of query sequences 
#TOT = 27143 #TODO derive from input
###

# read in genetic map, if supplied
# format for spruce map is LG\tcM\tcDNA
if genetMap:
    mapDat = pd.read_csv(genetMap, sep="\t", comment='#', low_memory=False, header=None, names=['LG', 'cM', 'cDNA'])
    # limit genetic map analysis to complete (i.e. single) cDNAs to improve confidence
    map_cDNA = set(mapDat.cDNA.tolist())
    uniqDatMap = uniqDat[uniqDat.qname.isin(map_cDNA)]
    uMap = uniqDatMap[uniqDatMap.tname.duplicated(keep=False)]

# setup counters
num_complete = 0
num_duplicated = 0
num_partial = 0
num_fragmented = 0
num_poor = 0

# apply check_complete to whole set
for rec in uniqDat.itertuples():
    res = check_aln(rec, 'assess')
    if res == 'complete':
        num_complete += 1
    elif res == 'partial':
        num_partial += 1
    elif res == 'poorly mapped':
        num_poor += 1

# apply check_frag to whole set
for qry in tlocDat.qname.unique():
    this_qry = tlocDat[tlocDat.qname == qry]
    frags = []
    for rec in this_qry.itertuples():
        frags.append(rec)
    
    res = check_frag(frags)
    if res == 'fragmented':
        num_fragmented += 1
    elif res == 'partial':
        num_partial += 1
    elif res == 'poorly mapped':
        num_poor += 1

# apply check_dupl to whole set
for qry in duplDat.qname.unique():
    this_qry = duplDat[duplDat.qname == qry]
    frags = []
    for rec in this_qry.itertuples():
        frags.append(rec)
    
    res = check_dupl(frags)    
    if res == 'complete':
        num_complete += 1
    elif res == 'duplicated':
        num_duplicated += 1
    elif res == 'partial':
        num_partial += 1
    elif res == 'poorly mapped':
        num_poor += 1

# calc percentages and report results
rate_complete = float(num_complete) / float(TOT)
rate_duplicated = float(num_duplicated) / float(TOT)
rate_partial = float(num_partial) / float(TOT)
rate_fragmented = float(num_fragmented) / float(TOT)
rate_poor = float(num_poor) / float(TOT)

pct_complete = round(100.0 * rate_complete, 2)
pct_duplicated = round(100.0 * rate_duplicated, 2)
pct_partial = round(100.0 * rate_partial, 2)
pct_fragmented = round(100 * rate_fragmented, 2)
pct_poor = round(100 * rate_poor, 2)

# count how many sequences are missing alignments
obs = count_aligned(uniqDat, duplDat, tlocDat)
num_miss = TOT - obs
rate_miss = float(num_miss) / float(TOT)
pct_miss = round(100.0 * rate_miss, 2)

# determine if the right number of sequences have a reported result
num_counted = sum([num_complete, num_duplicated, num_fragmented, num_partial, num_poor, num_miss])
rate_counted = float(num_counted) / float(TOT)
pct_counted = round(100 * rate_counted, 2)

# write to tsv
tsvout = "-".join([prefix, "results.tsv"])
with open(tsvout, "w") as outfile:
    header = "\t".join(["", "Complete", "Duplicated", "Fragmented", "Partial", "Poorly Mapped", "Missing", "Total cDNAs searched"])
    nums = "\t".join([str(x) for x in ["Number", num_complete, num_duplicated, num_fragmented, num_partial, num_poor, num_miss, num_counted]])
    pcts = "\t".join([str(x) for x in ["Percent", pct_complete, pct_duplicated, pct_fragmented, pct_partial, pct_poor, pct_miss, pct_counted]])

    print >> outfile, header
    print >> outfile, nums
    print >> outfile, pcts

jiraout = "-".join([prefix, "results.jira"])
with open(jiraout, "w") as outfile:
    header = "||".join(["", "Complete", "Duplicated", "Fragmented", "Partial", "Poorly Mapped", "Missing", "Total cDNAs searched", ""])
    nums = [num_complete, num_duplicated, num_fragmented, num_partial, num_poor, num_miss, num_counted]
    pcts = [pct_complete, pct_duplicated, pct_fragmented, pct_partial, pct_poor, pct_miss, pct_counted]
    res = "|" + "|".join([jira_formatter(x) for x in zip(nums, pcts)]) + "|"
    
    print >> outfile, header
    print >> outfile, res

# print to STDOUT
print "\n=== GNAVIGATOR RESULTS ==="
print "%s (%s%%) complete sequences" % (num_complete, pct_complete)
print "%s (%s%%) duplicated sequences" % (num_duplicated, pct_duplicated)
print "%s (%s%%) fragmented sequences" % (num_fragmented, pct_fragmented)
print "%s (%s%%) partial sequences" % (num_partial, pct_partial)
print "%s (%s%%) poorly mapped sequences" % (num_poor, pct_poor)
print "%s (%s%%) missing sequences" % (num_miss, pct_miss)
print "%s (%s%%) sequences were evaluated" % (num_counted, pct_counted)

# apply check_LG to whole uniq set
if genetMap:
    num_goodLG = 0 # same LG, right order
    num_WO_LG = 0 # same LG, wrong order
    num_diffLG = 0 # different LG

    for rec in uMap.tname.unique():
        thisScaf = uMap[uMap.tname.isin([rec])]
        res = check_LG(thisScaf, mapDat)
        if res == 'same LG, right order':
            num_goodLG += 1
        elif res == 'same LG, wrong order':
            num_WO_LG += 1
        elif res == 'different LG':
            num_diffLG += 1

# report genetic map results
if genetMap:
    num_scaff_toCheck = len(uMap.tname.unique())
    num_scaff_checked = num_goodLG + num_WO_LG + num_diffLG
    if num_scaff_toCheck == num_scaff_checked:
        rate_LGscaff = float(num_scaff_checked) / float(TOT)
        rate_goodLG = float(num_goodLG) / float(num_scaff_checked)
        rate_WO_LG = float(num_WO_LG) / float(num_scaff_checked)
        rate_diffLG = float(num_diffLG) / float(num_scaff_checked)

        pct_LGscaff = round(100.0 * rate_LGscaff, 2) 
        pct_goodLG = round(100.0 * rate_goodLG, 2)
        pct_WO_LG = round(100.0 * rate_WO_LG, 2)
        pct_diffLG = round(100.0 * rate_diffLG, 2)

    else:
        print 'Not all scaffolds to be checked against genetic map were successfully checked.'
        print 'Maybe something is wrong with the input data?'
        sys.exit(2)
    
    # write to tsv
    tsvout = "-".join([prefix, "genetic-map-results.tsv"])
    with open(tsvout, "w") as outfile:
        header = "\t".join(["", "Same LG, right order", "Same LG, wrong order", "Different LG", "Total scaffolds analyzed"])
        nums = "\t".join([str(x) for x in ["Number", num_goodLG, num_WO_LG, num_diffLG, num_scaff_checked]])
        pcts = "\t".join([str(x) for x in ["Percent", pct_goodLG, pct_WO_LG, pct_diffLG, pct_LGscaff]])

        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts

    jiraout = "-".join([prefix, "genetic-map-results.jira"])
    with open(jiraout, "w") as outfile:
        header = "||".join(["", "Same LG, right order", "Same LG, wrong order", "Different LG", "Total scaffolds analyzed", ""])
        nums = [num_goodLG, num_WO_LG, num_diffLG, num_scaff_checked]
        pcts = [pct_goodLG, pct_WO_LG, pct_diffLG, pct_LGscaff]
        res = "|" + "|".join([jira_formatter(x) for x in zip(nums, pcts)]) + "|"

        print >> outfile, header
        print >> outfile, res

    print "\n=== GENETIC MAP GNAVIGATOR RESULTS ==="
    print "%s (%s%%) scaffolds had 2+ complete cDNAs from the genetic map aligned to them." % (num_scaff_checked, pct_LGscaff)
    print "%s (%s%%) case(s) were from the same linkage group and in the expected order." % (num_goodLG, pct_goodLG)
    print "%s (%s%%) case(s) were from the same linkage group, but NOT in the expected order." % (num_WO_LG, pct_WO_LG)
    print "%s (%s%%) case(s) were from different linkage groups." % (num_diffLG, pct_diffLG)

### EOF ###
