"""
Modified from code written by brentp
Retrieved from https://www.biostars.org/p/710/ on May 5, 2015
"""

from itertools import groupby


def fasta_iter(fasta_name):
    """given a fasta file. yield tuples of header, sequence"""
    # ditch the boolean (x[0]) and just keep the header
    # or sequence since we know they alternate.
    faiter = (x[1] for x in groupby(fasta_name,
              lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

