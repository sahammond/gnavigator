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
            return (cDNA, scaf, 'Complete')
        elif pid >= ident_thold and pcov < cov_thold:
            if pcov >= 0.5:
                return (cDNA, scaf, 'Partial')
            else:
                return (cDNA, scaf, 'Poorly mapped')
        else:
            return (cDNA, scaf, 'Poorly mapped')
    
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
    scaf = []
        
    for aln in alns:
        this_aln = check_aln(aln, 'report', ident_thold, cov_thold)
        # check if if might qualify as complete
        check_comp = check_aln(aln, 'assess', ident_thold, cov_thold)
        if check_comp == 'Complete':
            scaf = this_aln[5]
            return (cDNA, scaf, 'Complete')
        goodb += this_aln[0]
        covb += this_aln[1]
        seg += this_aln[2]
        scaf.append(this_aln[5])
    else:
        qsize = this_aln[3]
        cDNA = this_aln[4]
    
    pid = goodb / seg
    pcov = covb / qsize
    
    scaf_rep = ";".join(scaf)
    if pid >= ident_thold and pcov >= cov_thold:
        return (cDNA, scaf_rep, 'Fragmented')
    elif pid >= ident_thold and pcov < cov_thold:
        if pcov > 0.5:
            return (cDNA, scaf_rep, 'Partial')
        else:
            return (cDNA, scaf_rep, 'Poorly mapped')
    else:
        return (cDNA, scaf_rep, 'Poorly mapped')


def check_dupl(alns, ident_thold, cov_thold):
    """check if a sequence has 2+ complete alignments"""
    # >1 alignment must be complete in order to be duplicated
    # one or more partial, fragmented, poorly mapped will not count
    cDNA = ''
    scaf = []
    results = {'Complete':[], 'Partial':[], 'Duplicated':[],
               'Poorly mapped':[]}
    for aln in alns:
        this_aln = check_aln(aln, 'assess', ident_thold, cov_thold)
        res = this_aln[-1]
        results[res].append(this_aln)
        scaf.append(this_aln[1])
    else:
        cDNA = this_aln[0]
    scaf_rep = ";".join(scaf)
    num_complete = len(results['Complete'])
    if num_complete == 1:
        best_scaf = scaf[0]
        return (cDNA, best_scaf, 'Complete')
    elif num_complete > 1:
        return (cDNA, scaf_rep, 'Duplicated')
    else:
        if len(results['Partial']) >= 1:
            return (cDNA, scaf_rep, 'Partial')
        else:
            return (cDNA, scaf_rep, 'Poorly mapped')
