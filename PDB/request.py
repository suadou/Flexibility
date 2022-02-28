"""
This module search in PDB and SwissProt database PDB files from X-ray difraction and AlphaFold prediction.
It returns the PDB file.
"""

from Bio.Blast import NCBIWWW
import prody
import requests
import re


def BLAST_PDB(query):
    """
    Using a query class as input it returns a text file with the matches from BLASTp search in PDB database,
    """
    blast = prody.blastPDB(query.sequence)
    best_hit = blast.getBest()
    return best_hit


def dowload_PDB(query, pdb_id):
    """
    Download PDB files from PDB.
    """
    url = "https://files.rcsb.org/download/" + \
        pdb_id['pdb_id'].upper() + ".pdb"
    blast = requests.get(url)
    fo = open(query.identifier+'.pdb', 'w')
    for line in blast:
        fo.write(line)
    fo.close()


def BLAST_sp(query):
    blast_XML = NCBIWWW.qblast("blastp", database="swissprot", sequence=query.sequence,
                               alignments=0, hitlist_size=50, expect=1e-5, format_type="XML")
    return blast_XML


def parse_XML(blast_XML):
    swissprot_id = ''
    hsp_alignlen = ''
    hsp_identity = ''
    for line in blast_XML:
        if re.search('<Hit_id>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hit_id>')
            line = line.rstrip('</')
            line = line.strip('sp|')
            line = line.rstrip('|')
            swissprot_id = line.split('.')[0]
        if re.search('<Hsp_identity>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hsp_identity>')
            hsp_identity = line.rstrip('</')
        if re.search('<Hsp_align-len>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hsp_align-len>')
            hsp_alignlen = line.rstrip('</')
        if swissprot_id and hsp_alignlen and hsp_identity:
            tuple([swissprot_id, int(hsp_identity) / int(hsp_alignlen)])
            swissprot_id = ''
            hsp_alignlen = ''
            hsp_identity = ''


def download_AlphaFold(query, uniprot_id):
    url = "https://alphafold.ebi.ac.uk/files/AF-" + uniprot_id + "-F1-model_v2.pdb"
    pdb = requests.get(url)
    fo = open(query.identifier+'_AlphaFold.pdb', 'w')
    for line in pdb:
        fo.write(line)
    fo.close()
