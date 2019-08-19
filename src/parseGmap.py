"""
Created on 22 May 2018

@author: Lauren Coombe

Given gnavigator output and an input genetic map file, output a file of the following format:
LG\tcMg\tcDNA\tscafID+12030304
Adding scaffold ID (clone cDNA best aligns to), orientation, scaffold size
"""

#from __future__ import print_function
#import argparse


#Storing information about a scaffold
class scaffoldInfo:
    
    def __init__(self, id, orient=None, seqID=float("-inf"), cov=float("-inf"), size=float("-inf")):
        self.id = id
        self.orient = orient
        self.seqID = seqID
        self.cov = cov
        self.size = size

#Read in the gnavigator cDNA full table TSV results, return dictionary of cDNA --> [scaffoldDicts]
def readGnavResults(gnavOut):
    cDNAs = {}
    with open(gnavOut, 'r') as gnav:
        header = gnav.readline()
        for line in gnav:
            line = line.strip().split("\t")
            (cDNA, scaf, lg) = (line[0], line[2], line[3])
            if lg == "unknownLG" or scaf == "NA":
                continue
            if cDNA not in cDNAs:
                cDNAs[cDNA] = {}
            scafInfo = scaffoldInfo(scaf)
            cDNAs[cDNA][scaf] = scafInfo
#    print 'cDNAs ' + str(len(cDNAs)) ###
    return cDNAs

#Reads through gmap PSL alignments, track the seqID, cov, orient information for each needed scaffold
def readAligns(cDNAs, aligns_file):
#    print aligns_file ###
    with open(aligns_file, 'r') as aligns:
        for align in aligns:
            align = align.strip().split("\t")
            (matches, mismatches, qinserts, strand) = (int(align[0]), int(align[1]), int(align[4]), align[8])
            (qname, qstart, qend, qsize) = (align[9], int(align[11]), int(align[12]), int(align[10]))
            (tname, tsize) = (align[13], int(align[14]))
            cov = (matches + mismatches - qinserts)/float(qsize)
            seqID = float(matches)/(qend - qstart)

            #Only track information if it is one of the target scaffolds from gnavigator
            if qname in cDNAs and tname in cDNAs[qname]:
                if cDNAs[qname][tname].cov is None or cDNAs[qname][tname].seqID is None\
                        or cov > cDNAs[qname][tname].cov \
                        or (cDNAs[qname][tname].cov == cov and seqID > cDNAs[qname][tname].seqID):
                    cDNAs[qname][tname] = scaffoldInfo(tname, strand, seqID, cov, tsize)

#Given a dictionary of scaffolds mapping to a cDNA, return the scaffoldInfo object of the 'best hit'
def getBestHit(scaffolds):
    curBest = None
    curBestCov = float("-inf")
    curBestSeqID = float("-inf")
#    print scaffolds ###
    for scaf in scaffolds:
        if scaffolds[scaf].cov > curBestCov \
            or (scaffolds[scaf].cov == curBestCov and scaffolds[scaf].seqID > curBestSeqID)\
            or (scaffolds[scaf].cov == curBestCov and scaffolds[scaf].seqID == curBestSeqID and scaffolds[scaf].size > curBest.size):
            curBest = scaffolds[scaf]
            curBestCov = scaffolds[scaf].cov
            curBestSeqID = scaffolds[scaf].seqID
#    print "Current best scaf is %s" % curBest.id ###
    return curBest


#Read through the ordered genetic map, adding 4th column with extra information
def outputGM(cDNAs, GM_file):
    with open(GM_file, 'r') as GM:
        for line in GM:
            line = line.strip().split("\t")
            (lg, cMg, cDNA) = (line[0], line[1], line[2])
            if cDNA in cDNAs:
                # Output the best scaffold hit
#                print cDNA ###
                scaffolds = cDNAs[cDNA]
                scaf = getBestHit(scaffolds)
                scafStr = "%s%s%s" % (scaf.id, scaf.orient, scaf.size)
#               print("%s\t%s\t%s\t%s" % (lg, cMg, cDNA, scafStr))
                outbuff = "%s\t%s\t%s\t%s" % (lg, cMg, cDNA, scafStr)
                yield outbuff
            else:
                #print("%s\t%s\t%s\tmissing" % (lg, cMg, cDNA))
                outbuff = "%s\t%s\t%s\tmissing" % (lg, cMg, cDNA)
                yield outbuff

def wrapper(gnavOut, uniq, mult, tloc, GM, prefix):
#   parser = argparse.ArgumentParser(description="Get information about genetic map coverage scale")
#   parser.add_argument("gnavOut", help="Gnavigator full cDNA results TSV table")
#   parser.add_argument("GM", help="Genetic map")
#   parser.add_argument("-u", help="Gmap .uniq alignment file", required=True)
#   parser.add_argument("-m", help="Gmap .mult alignment file", required=True)
#   parser.add_argument("-t", help="Gmap .transloc alignment file", required=True)

#   args = parser.parse_args()

    cDNAs = readGnavResults(gnavOut)
    readAligns(cDNAs, uniq)
    readAligns(cDNAs, mult)
    readAligns(cDNAs, tloc)

    expanded_gm = '-'.join([prefix, 'expanded-genetic-map.tsv'])
    with open(expanded_gm, 'w') as outfile:
        for line in outputGM(cDNAs, GM):
            print(line, file=outfile)
#outputGM(cDNAs, GM)


#if __name__ == '__main__':
#   main()
