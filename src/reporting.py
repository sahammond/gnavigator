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
    missingL = [(x, 'NA', 'Missing') for x in missing]
    
    return (missingL, tot_cDNA)
