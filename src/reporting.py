"""
Created on 23 May 2018

@author: S. Austin Hammond

Functions to assess overall results.
"""

import sys
import os

# ensure gnavigator's src directory is in PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/gnavigator/src')

from src.fasta_iter import fasta_iter
import src.utilities as util

def find_missing(fasta, res_dict):
    # fasta is cDNA fasta
    # res_dict is results dictionary like cDNA_res
    cDNA_set = set()
    with open(fasta, 'r') as infile:
        for rec in fasta_iter(infile):
            seqid = rec[0]
            seqid_only = seqid.split(" ")[0]
            cDNA_set.add(seqid_only)
    tot_cDNA = len(cDNA_set)
    detected = set()
    for key, value in res_dict.items():
        resID = set([x[0] for x in value])
        detected = detected.union(resID)
    missing = cDNA_set.difference(detected)
    missingL = [(x, 'NA', 'Missing', 'NA', 'NA') for x in missing]
    
    return (missingL, tot_cDNA)


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

    sum_complete = num_complete + num_duplicated
    pct_sum_complete = pct_complete + pct_duplicated

    # report if the right number of sequences have a result
    num_counted = sum([num_complete, num_duplicated, num_fragmented,
                       num_partial, num_poor, num_missing])
    rate_counted = float(num_counted) / float(TOT)
    pct_counted = round(100 * rate_counted, 2)

    # write to tsv
    tsvout = '-'.join([prefix, 'results.tsv'])
    header_txt = ['', 'Complete', 'Complete, single copy', 'Complete, multiple copies',
                  'Fragmented', 'Partial', 'Poorly Mapped', 'Missing', 'Total cDNAs searched']
    with open(tsvout, 'w') as outfile:
        header = '\t'.join(header_txt)
        nums = '\t'.join([str(x) for x in ['Number', sum_complete, num_complete, num_duplicated,
                                           num_fragmented, num_partial, num_poor,
                                           num_missing, num_counted]])
        pcts = '\t'.join([str(x) for x in ['Percent', pct_sum_complete, pct_complete,
                                           pct_duplicated, pct_fragmented, pct_partial, pct_poor,
                                           pct_missing, pct_counted]])
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts
    jiraout = '-'.join([prefix, 'results.jira'])
    with open(jiraout, 'w') as outfile:
        header_txt.extend(['']) # add final empty string for jira table separator
        header = '||'.join(header_txt)
        nums = [sum_complete, num_complete, num_duplicated, num_fragmented, num_partial, num_poor,
                num_missing, num_counted]
        pcts = [pct_sum_complete, pct_complete, pct_duplicated, pct_fragmented, pct_partial,
                pct_poor, pct_missing, pct_counted]
        res = '|' + '|'.join([util.jira_formatter(x) for x in zip(nums, pcts)]) + '|'
        
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, res

    # print to STDOUT
    print '\n=== GNAVIGATOR cDNA RESULTS ==='
    print '%s (%s%%) complete sequences' % (sum_complete, pct_sum_complete)
    print '  %s (%s%%) complete, single copy sequences' % (num_complete, pct_complete)
    print '  %s (%s%%) complete, multiple copy sequences' % (num_duplicated, pct_duplicated)
    print '%s (%s%%) fragmented sequences' % (num_fragmented, pct_fragmented)
    print '%s (%s%%) partial sequences' % (num_partial, pct_partial)
    print '%s (%s%%) poorly mapped sequences' % (num_poor, pct_poor)
    print '%s (%s%%) missing sequences' % (num_missing, pct_missing)
    print '%s (%s%%) sequences were evaluated' % (num_counted, pct_counted)
    print util.report_time()


def output_cDNA(prefix, cDNA_results, gm=''):
    # write out cDNA:scaffold mappings
    full_out = '-'.join([prefix, 'full-cDNA-results-table.tsv'])
    with open(full_out, 'w') as outfile:
        if len(gm) > 0:
            header = '\t'.join(['# cDNA ID', 'Status', 'Scaffold',
                                'Linkage group'])
            print >> outfile, util.report_cmd()
            print >> outfile, header
            for status, result in cDNA_results.items():
                for res in result:
                    for t in util.table_formatter_wGM(res, gm):
                        print >> outfile, t
        else:
            header = '\t'.join(['# cDNA ID', 'Status', 'Scaffold'])
            print >> outfile, util.report_cmd()
            print >> outfile, header
            for status, result in cDNA_results.items():
                for res in result:
                    for t in util.table_formatter(res):
                        print >> outfile, t


# write out cDNA:scaffold mappings
def output_gm(prefix, gm_res):
    header = '\t'.join(['# Scaffold', 'cDNA IDs', 'Status', 'Linkage group(s)'])
    full_out = '-'.join([prefix, 'full-genetic-map-results-table.tsv'])
    with open(full_out, 'w') as outfile:
        print >> outfile, util.report_cmd()
        print >> outfile, header
        for status, result in gm_res.items(): ### what do I do with status here? Or is it just discarded because there should only be complete, single cDNAs in gm_res
            for res in result:
                for t in util.LG_table_formatter(res):
                    print >> outfile, t


def report_gm(uniqDatMap_select, gm_results, gm_cdna_statuses, prefix):
    """report summary of genetic map results"""
    num_scaff_toCheck = len(uniqDatMap_select.tname.unique())
    num_solo = gm_cdna_statuses['Solo']
    #tot = num_scaff_toCheck + gm_cdna_statuses['Solo']
    num_goodLG = len(gm_results['goodLG'])
    num_WO_LG = len(gm_results['WO_LG'])
    num_diffLG = len(gm_results['diffLG'])
    num_undet = len(gm_results['undet'])
    num_scaff_checked = num_goodLG + num_WO_LG + num_diffLG + num_undet + num_solo
    num_2plus_scaff = num_goodLG + num_WO_LG + num_diffLG + num_undet

    tot = num_scaff_checked # for legacy reasons
    if num_scaff_toCheck == num_scaff_checked:
        rate_solo = float(num_solo) / float(tot)
        rate_LGscaff = float(num_2plus_scaff) / float(tot)
        rate_goodLG = float(num_goodLG) / float(tot)
        rate_WO_LG = float(num_WO_LG) / float(tot)
        rate_diffLG = float(num_diffLG) / float(tot)
        rate_undet = float(num_undet) / float(tot)
        pct_solo = round(100.0 * rate_solo, 2)
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
    tsvout = '-'.join([prefix, 'genetic-map-scaffold-results.tsv'])
    with open(tsvout, 'w') as outfile:
        header_txt = ['', '1 complete cDNA', '2+ complete cDNAs', 'Same LG, expected order',
                      'Same LG, unexpected order','Different LG', 'Same LG, undetermined order']
        header = '\t'.join(header_txt)
        nums = '\t'.join([str(x) for x in ['Number', num_solo, num_2plus_scaff, num_goodLG,
                                           num_WO_LG, num_diffLG, num_undet]])
        pcts = '\t'.join([str(x) for x in ['Percent', pct_solo, pct_LGscaff, pct_goodLG,
                                           pct_WO_LG, pct_diffLG, pct_undet]])
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts
    jiraout = '-'.join([prefix, 'genetic-map-scaffold-results.jira'])
    with open(jiraout, 'w') as outfile:
        header_txt.extend(['']) # add final empty string for jira table separator
        header = '||'.join(header_txt)
        nums = [num_solo, num_2plus_scaff, num_goodLG, num_WO_LG, num_diffLG, num_undet]
        pcts = [pct_solo, pct_LGscaff, pct_goodLG, pct_WO_LG, pct_diffLG, pct_undet]
        res = '|' + '|'.join([util.jira_formatter(x) for x in zip(nums, pcts)]) + '|'
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, res

    # print to STDOUT
    print '\n=== GNAVIGATOR GENETIC MAP RESULTS per SCAFFOLD with 1+ COMPLETE cDNA ==='
    print ' '.join(['%s (%s%%) had exactly 1 complete cDNA from the genetic map',
                    'aligned to them.']) % (num_solo, pct_solo)
    print ' '.join(['%s (%s%%) had 2+ complete cDNAs from the genetic map',
                    'aligned to them.']) % (num_2plus_scaff, pct_LGscaff)
    print ' '.join(['  %s (%s%%) were from the same linkage group and in the',
                    'expected order.']) % (num_goodLG, pct_goodLG)
    print ' '.join(['  %s (%s%%) were from the same linkage group, but NOT in',
                    'the expected order.']) % (num_WO_LG, pct_WO_LG)
    print '  %s (%s%%) were from different linkage groups.' % (num_diffLG, pct_diffLG)
    print ' '.join(['  %s (%s%%) were from the same linkage group but their order',
                    'could not be determined.']) % (num_undet, pct_undet)
    print util.report_time()


def measure_scaf(fasta):
    """Create dict of scaffold ID and size"""
    assembly = {}
    with open(fasta, 'r') as infile:
        for rec in fasta_iter(infile):
            seqid, sqn = rec
            seqid_only = seqid.split(" ")[0]
            slen = len(sqn)
            assembly[seqid_only] = slen

    return assembly
    

def report_gm_cDNA(gm_results, uniqMapDat_select, prefix):
    """Report genetic map results as a proportion of eligible cDNAs"""

    # input cDNA_results is dict with keys == cDNA status
    # see src.classify.check_aln
    gm_cdna_res = {}
    statuses = {'Solo':0, 'Same LG, expected order':0, 'Same LG, unexpected order':0,
                'Different LG':0, 'Same LG, order undetermined':0}
    tot = 0

    for res in gm_results.items():
        for pres in res[1]:
            scaf, cdnaS, stat, lg = pres
            cdnaL = cdnaS.split(';')
            for cdna in cdnaL:
                gm_cdna_res[cdna] = [scaf, stat]
                statuses[stat] += 1
                tot += 1

    miss = 'Solo'
    for rec in uniqMapDat_select.iterrows():
        cdna = rec[1][9]
        scaf = rec[1][13]

        if cdna not in gm_cdna_res:
            gm_cdna_res[cdna] = [scaf, miss]
            statuses[miss] += 1
            tot += 1

    num_solo = statuses['Solo']
    num_good = statuses['Same LG, expected order']
    num_wo = statuses['Same LG, unexpected order']
    num_diff = statuses['Different LG']
    num_undet = statuses['Same LG, order undetermined']

    rate_solo = float(statuses['Solo']) / float(tot)
    rate_good = float(statuses['Same LG, expected order']) / float(tot)
    rate_wo = float(statuses['Same LG, unexpected order']) / float(tot)
    rate_diff = float(statuses['Different LG']) / float(tot)
    rate_undet = float(statuses['Same LG, order undetermined']) / float(tot)

    pct_solo = round(100.0 * rate_solo, 2)
    pct_good = round(100.0 * rate_good, 2)
    pct_wo = round(100.0 * rate_wo, 2)
    pct_diff = round(100.0 * rate_diff, 2)
    pct_undet = round(100.0 * rate_undet, 2)

    # write to tsv
    tsvout = '-'.join([prefix, 'genetic-map-cDNA-results.tsv'])
    with open(tsvout, 'w') as outfile:
        header_txt = ['', 'Solo (unevaluated)', 'Same LG, expected order',
                      'Same LG, unexpected order', 'Different LG', 'Same LG, order undetermined']
        header = '\t'.join(header_txt)
        nums = '\t'.join([str(x) for x in ['Number', num_solo, num_good, num_wo, num_diff,
                                           num_undet]])
        pcts = '\t'.join([str(x) for x in ['Percent', pct_solo, pct_good, pct_wo, pct_diff,
                                           pct_undet]])
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, nums
        print >> outfile, pcts
    jiraout = '-'.join([prefix, 'genetic-map-cDNA-results.jira'])
    with open(jiraout, 'w') as outfile:
        header_txt.extend(['']) # add final empty string for jira table separator
        header = '||'.join(header_txt)
        nums = [num_solo, num_good, num_wo, num_diff, num_undet]
        pcts = [pct_solo, pct_good, pct_wo, pct_diff, pct_undet]
        res = '|' + '|'.join([util.jira_formatter(x) for x in zip(nums, pcts)]) + '|'
        print >> outfile, util.report_cmd()
        print >> outfile, header
        print >> outfile, res

    print '\n=== GNAVIGATOR GENETIC MAP RESULTS per COMPLETE cDNA ==='
    print ' '.join(['%s (%s%%) were found alone on a scaffold and were not evaluated',
                    'further.']) % (num_solo, pct_solo)
    print ' '.join(['%s (%s%%) were from the same linkage group and were found on scaffolds in',
                    'the expected order']) % (num_good, pct_good)
    print ' '.join(['%s (%s%%) were from the same linkage group but were found on scaffolds in an',
                    'unexpected order.']) % (num_wo, pct_wo)
    print ' '.join(['%s (%s%%) were from different linkage groups but were found on the same',
                    'scaffold.']) % (num_diff, pct_diff)
    print ' '.join(['%s (%s%%) were from the same linkage group but their order could not be',
                    'determined.']) % (num_undet, pct_undet)
    print util.report_time()

    # return statuses so certain stats can be used elsewhere
    return statuses
