"""
Created on 23 May 2018

@author: S. Austin Hammond

Helper and formatting functions.
"""

from time import localtime, strftime


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
