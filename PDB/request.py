"""
request module search in PDB and SwissProt database using BLASTp looking for PDB
files from X-ray difraction and AlphaFold prediction with the best percentage
of identity.
"""

from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
import requests
import re


def BLAST_commandline(query, database):
    blast_XML = NcbiblastpCommandline(cmd='blastp', query=query.sequence,
                                      db=database, alignments=0, hitlist_size=50, expect=1e-5, format_type="XML")
    return blast_XML


def BLAST_PDB(query):
    """
    Using a Query class as input it returns a XML BLASTp results from PDB
    database.
    """
    blast_XML = NCBIWWW.qblast("blastp", database="pdb", sequence=query.sequence,
                               alignments=0, hitlist_size=50, expect=1e-5, format_type="XML")
    return blast_XML


def download_PDB(query, pdb_id, path):
    """
    Download PDB files from PDB. It returns a pdb file ([identifier].pdb) in the
    path directory. The arguments are the following:
        query - It is a query class wich contains sequence and identifier.
        pdb_id - It is the PDB identifier.
        path - Directory to download the PDB file.
    """
    url = "https://files.rcsb.org/download/" + \
        pdb_id + ".pdb"
    if check_url(url):
        blast = requests.get(url)
        fo = open(path + query.identifier+'.pdb', 'w')
        for line in blast:
            fo.write(line.decode("utf-8"))
        fo.close()
    else:
        return False


def BLAST_sp(query):
    """
    Using a Query class as input it returns a XML BLASTp results from SwissProt
    database.
    """
    blast_XML = NCBIWWW.qblast("blastp", database="swissprot", sequence=query.sequence,
                               alignments=0, hitlist_size=50, expect=1e-5, format_type="XML")
    return blast_XML


def download_AlphaFold(query, swissprot_id, path):
    """
    Download PDB files from AlphaFold database. It returns a pdb file
    ([identifier]_AlphaFold.pdb) in the path directory. The arguments are the
    following:
        query - It is a query class wich contains sequence and identifier.
        swissprot_id - It is the SwissProt identifier.
        path - Directory to download the PDB file.
    """
    url = "https://alphafold.ebi.ac.uk/files/AF-" + swissprot_id + "-F1-model_v2.pdb"
    if check_url(url):
        pdb = requests.get(url)
        fo = open(path + query.identifier+'_AlphaFold.pdb', 'w')
        for line in pdb:
            fo.write(line.decode("utf-8"))
        fo.close()
    else:
        return False


def parse_XML(blast_XML):
    """
    parse_XML is a generator function which use a XML BLASTp results from NCBI
    as input to return a list containing PDB/SwissProt id, perecentage of identity
    the start and end position, and the gapas from each match.
    """
    swissprot_id = ''
    hsp_alignlen = ''
    hsp_identity = ''
    start = ''
    end = ''
    gaps = []
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
        if re.search('<Hsp_hit-from>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hsp_hit-from>')
            start = line.rstrip('</')
        if re.search('<Hsp_hit-to>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hsp_hit-to>')
            end = line.rstrip('</')
        if re.search('<Hsp_hseq>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = line.strip('<Hsp_hseq>')
            line = line.rstrip('</')
            for i in range(0, len(line)):
                if line[i] == '-':
                    gaps.append(i+1)
        if swissprot_id and hsp_alignlen and hsp_identity and start and end:
            identity = int(hsp_identity) / int(hsp_alignlen)
            yield (str(swissprot_id), identity, start, end, gaps)
            swissprot_id = ''
            hsp_alignlen = ''
            hsp_identity = ''
            identity = ''
            start = ''
            end = ''
            gaps = []


def check_url(url):
    """ Check if a url exist or not. It returns a boolean"""
    request_response = requests.head(url)
    status_code = request_response.status_code
    website_is_up = status_code == 200
    return website_is_up


class PDB(object):
    """
    PDB object stores the Query identifier of the record, the PDB identifier,
    the PDB file path downloaded, the percentage of identity of the match, the
    position of start and end of the match and the gap position.
    """

    def __init__(self, query, identifier, path, identity, start, end, gaps):
        self.query = query
        self.identifier = identifier
        self.path = path
        self.identity = identity
        self.start = start
        self.end = end
        self.gaps = gaps
