"""
Created on 23 May 2018

@author: S. Austin Hammond

Function to assess features from a genetic map.
"""

from itertools import permutations


def check_LG(query, genetic_map):
    """check if cDNAs are from the same LG"""
    # query is a pandas array of a single scaffold's alignments
    # genetic_map is a pandas array of LG\tcM\tcDNA
    # assume for now that cDNA is in genetic_map$cDNA
    # returns 'same LG, right order', 'same LG, wrong order',
    ## 'different LG', or 'order undetermined'
    refs = query.qname.tolist()
    thisMap = genetic_map[genetic_map.cDNA.isin(refs)]
    theseLGlist = thisMap.LG.unique().tolist()
    theseLG = ";".join(theseLGlist)
    numLG = len(theseLGlist)
    scaf = str(query.tname.unique()[0])
    fwdA = query.sort_values(['tstart'])
    fwdL = fwdA.qname.tolist()
    revA = query.sort_values(['tstart'], ascending=False)
    revL = revA.qname.tolist()
    mapL = thisMap.cDNA.tolist()
    cDNA_names = ";".join(fwdL)

    if len(thisMap) == 2:
        if numLG == 1:
            return (scaf, cDNA_names, 'Same LG, right order', theseLG)
        else:
            return (scaf, cDNA_names, 'Different LG', theseLG)
    else:
        if numLG == 1:
            # compare both forward and reverse orders
            if mapL == fwdL or mapL == revL:
                return (scaf, cDNA_names, 'Same LG, right order', theseLG)
            # some GM features have the same position in cM
            # shuffle such features and check each permutation
            # currently only handles one block of same-cM features
            ## I haven't seen more than that on a single scaffold yet
            else:
                dup_cm = thisMap[thisMap.duplicated('cM', False)]
                num_dup = len(dup_cm.cM.unique())
                dup_cDNAs = dup_cm.cDNA.unique().tolist()
                perm = permutations(dup_cDNAs)

                # if all features are at same cM, no way to know if order is correct
                # give it the benefit of the doubt for now
                # in future, may change to 'Same LG, order undetermined'
                if len(dup_cm) == len(fwdL):
                    return (scaf, cDNA_names, 'Same LG, right order', theseLG)
                elif num_dup > 1:
                    # 2+ blocks of same-cM features not handled well
                    return (scaf, cDNA_names, 'Same LG, order undetermined', theseLG)
                else:
                    ord_cDNA = thisMap.sort_values('cM').cDNA.tolist()
                    ord_cDNA_base = [x for x in ord_cDNA if x not in dup_cDNAs]
                    # determine position in ord_cDNA_base to insert
                    ## the permuted same-cM cDNAs at
                    ins_pos = 0
                    for res in enumerate(ord_cDNA):
                        num = res[0]
                        nam = res[1]
                        if nam in dup_cDNAs:
                            ins_pos = num
                            break
                    # generate the various orders of cDNAs to check
                    permList = []
                    for permut in perm:
                        ord_cDNA_toCheck = ord_cDNA_base[:ins_pos] # preceding
                        ord_cDNA_toCheck.extend(permut) # permuted
                        ord_cDNA_toCheck.extend(ord_cDNA_base[ins_pos:]) # proceeding
                        permList.append(ord_cDNA_toCheck)                 
                    # check the orders
                    for rec in permList:
                        if fwdL == rec or revL == rec:
                            return (scaf, ";".join(rec), 'Same LG, right order', theseLG)
                    # if none of the permutations match, return 'wrong order'
                    else:
                        return (scaf, cDNA_names, 'Same LG, wrong order', theseLG)
        else:
            return (scaf, cDNA_names, 'Different LG', theseLG)
