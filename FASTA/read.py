"""
read module read a FASTA or multi-FASTA protein file and returns a Query object
for each sequence introduced.
"""
from Bio import SeqIO


class Query(object):
    """
    Query object store the identifier and the sequence of a FASTA protein
    sequence.
    """

    def __init__(self, sequence, identifier):
        self.sequence = sequence
        self.identifier = identifier


def FASTA_parser(input_file):
    """
    FASTA_parser use a FASTA or multi-FASTA protein file path and returns a list
    of Query objects of each record. 
    """
    fasta_sequence = SeqIO.parse(open(input_file), 'fasta')
    query_list = []
    for record in fasta_sequence:
        query_list.append(Query(record.seq, record.id))
    return query_list
