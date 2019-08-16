"""
Created on 23 May 2018

@author: S. Austin Hammond

These functions classify the alignments.

"""


def check_aln(aln, mode, ident_thold, cov_thold):
    """check the alignment of a cDNA sequence"""
    matches = float(aln.matches)
    mismatches = float(aln.mismatches)
    qinserts = float(aln.qbaseinsert)
    qsize = float(aln.qsize)
    qstart = float(aln.qstart)
    qend = float(aln.qend)
    scaf = str(aln.tname)
    cDNA = str(aln.qname)
    
    seg = qend - qstart
    pid = matches / seg
    # penalize insertions b/c reflects correctness of assembly
    pcov = (matches + mismatches - qinserts) / qsize
    
    if mode == 'assess':
        if pid >= ident_thold and pcov >= cov_thold:
            return (cDNA, scaf, 'Complete', pid, pcov)
        elif pid >= ident_thold and pcov < cov_thold:
            if pcov >= 0.5:
                return (cDNA, scaf, 'Partial', pid, pcov)
            else:
                return (cDNA, scaf, 'Poorly mapped', pid, pcov)
        else:
            return (cDNA, scaf, 'Poorly mapped', pid, pcov)
    
    elif mode == 'report':
        goodb = matches
        # penalize insertions b/c reflects correctness of assembly
        covb = matches + mismatches - qinserts
        
    return (goodb, covb, seg, qsize, cDNA, scaf)


def check_frag(alns, ident_thold, cov_thold):
    """check if 'chimeric' alignment is good or not"""
    # take in a series of alignments and evaluate them together
    # TODO extend to assess non-chimeric alignments
    ## some may be fragmented along >2 scaffolds
    ## with some missing bits in the middle
    goodb = 0
    covb = 0
    seg = 0
    qsize = 0
    cDNA = ''
    scafL = []
        
    for aln in alns:
        # check if it might qualify as complete
        cDNA, scaf, status, pid, pcov = check_aln(aln, 'assess', ident_thold, cov_thold)
        if status == 'Complete':
            return (cDNA, scaf, 'Complete', pid, pcov)
        this_aln = check_aln(aln, 'report', ident_thold, cov_thold)
        goodb += this_aln[0]
        covb += this_aln[1]
        seg += this_aln[2]
        scafL.append(this_aln[5])
    else:
        qsize = this_aln[3]
        cDNA = this_aln[4]
    
    pid = goodb / seg
    pcov = covb / qsize
    
    scaf_rep = ";".join(scafL)
    if pid >= ident_thold and pcov >= cov_thold:
        return (cDNA, scaf_rep, 'Fragmented', pid, pcov)
    elif pid >= ident_thold and pcov < cov_thold:
        if pcov > 0.5:
            return (cDNA, scaf_rep, 'Partial', pid, pcov)
        else:
            return (cDNA, scaf_rep, 'Poorly mapped', pid, pcov)
    else:
        return (cDNA, scaf_rep, 'Poorly mapped', pid, pcov)


def check_dupl(alns, ident_thold, cov_thold):
    """check if a sequence has 2+ complete alignments"""
    # >1 alignment must be complete in order to be duplicated
    # one or more partial, fragmented, poorly mapped will not count
    cDNA = ''
    scafL = []
    results = {'Complete':[], 'Partial':[], 'Poorly mapped':[]}
    for aln in alns:
        this_aln = check_aln(aln, 'assess', ident_thold, cov_thold)
        cDNA, scaf, status, pid, pcov = this_aln
        results[status].append(this_aln)
        scafL.append(scaf)
#    else:
#        cDNA = this_aln[0]
    scaf_rep = ";".join(scafL)
    num_complete = len(results['Complete'])
    if num_complete == 1:
        best_scaf = scafL[0] # best alignment appears first
        return (cDNA, best_scaf, 'Complete', pid, pcov)
    elif num_complete > 1:
        # get the best (first) one
        cdna, scaf, stat, pid, pcov = results['Complete'][0]
        return (cDNA, scaf_rep, 'Duplicated', pid, pcov)
    else:
        if len(results['Partial']) >= 1:
            # get the best (first) one
            cdna, scaf, stat, pid, pcov = results['Partial'][0]
            return (cDNA, scaf_rep, 'Partial', pid, pcov)
        else:
            # get the best (first) one
            cdna, scaf, stat, pid, pcov = results['Poorly mapped'][0]
            return (cDNA, scaf_rep, 'Poorly mapped', pid, pcov)


def assess(checkU, checkD, checkM, uniqDat, tlocDat, duplDat, cDNA_results,
           ident_thold, cov_thold):
    # apply check_complete to whole set
    if checkU:
        for rec in uniqDat.itertuples():
            res = check_aln(rec, 'assess', ident_thold, cov_thold)
            cDNA_results[res[2]].append(res) # append results tuple

    # apply check_frag to whole set
    if checkD:
        for qry in tlocDat.qname.unique():
            this_qry = tlocDat[tlocDat.qname == qry]
            frags = []
            for rec in this_qry.itertuples():
                frags.append(rec)
            res = check_frag(frags, ident_thold, cov_thold)
            cDNA_results[res[2]].append(res)

    # apply check_dupl to whole set
    if checkM:
        for qry in duplDat.qname.unique():
            this_qry = duplDat[duplDat.qname == qry]
            frags = []
            for rec in this_qry.itertuples():
                frags.append(rec)
            res = check_dupl(frags, ident_thold, cov_thold)
            cDNA_results[res[2]].append(res)

    return cDNA_results


def selex(cDNA_results):
    # lazy function to select cDNAs that are complete
    out = set()
    for res in cDNA_results['Complete']:
        this_cDNA = res[0]
        out.add(this_cDNA)

    return out
