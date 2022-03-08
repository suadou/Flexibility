import argparse
from PDB import request
from FASTA import read
import threading
import configparser


def AmphiProt(query, local, database=None):
    if local == False:
        BLAST = request.BLAST_sp(query)
    else:
        BLAST = request.BLAST_commandline(query, database)
    for seq in request.parse_XML(BLAST):
        if request.download_AlphaFold(query, seq[0], './files/PDB/') == False:
            continue
        else:
            AlphaFold_list.append(request.PDB(
                query.identifier, seq[0], './files/PDB/' + query.identifier + '_AlphaFold.pdb', seq[1], seq[2], seq[3], seq[4]))
            break


def PDB(query, local, database=None):
    if local == False:
        BLAST = request.BLAST_PDB(query)
    else:
        BLAST = request.BLAST_commandline(query, database)
    for seq in request.parse_XML(BLAST):
        if request.download_PDB(query, seq[0], './files/PDB/') == False:
            continue
        else:
            PDB_list.append(request.PDB(
                query.identifier, seq[0], './files/PDB/' + query.identifier + '.pdb', seq[1], seq[2], seq[3], seq[4]))
            break


parser = argparse.ArgumentParser(
        description="This program develops a flexibility score for proteins. It retreats a parseable text file and a graphycal representation of flexibility scores.")
parser.add_argument('-i', '-in', '--input', dest='input_file', action='store',
                    help='FASTA or Multi-FASTA protein sequence input file', required=True, default=None)
parser.add_argument('-o', '-out', '--output', dest='output_file', action='store',
                    help='Protein flexibility output file. It will give 2 files ([].txt - Parseable text file) and ([].png - Graphycal representation of flexibility scores)', required=False, default='flexibility_output')
parser.add_argument('--output_txt', dest='output_txt_file', action='store',
                    help='Protein flexibility parseable text file output.', required=False, default='flexibility_txt_output')
parser.add_argument('--output_png', dest='output_png_file', action='store',
                    help='Protein flexibility graphycal representation file.', required=False, default='flexibility_png_output')
options = parser.parse_args()

config = configparser.ConfigParser()

print(f"Reading {options.input_file} file.")
fasta_list = read.FASTA_parser(options.input_file)
print(f"{len(fasta_list)} record/s was/were read.")

# Threading for API BLASTp search
if config['blast']['local'] == False:
    threads = len(fasta_list)
    jobs = []
    PDB_list = []
    AlphaFold_list = []
    for i in range(0, threads):
        thread = threading.Thread(target=AmphiProt(fasta_list[i], False))
        jobs.append(thread)
        thread = threading.Thread(target=PDB(fasta_list[i], False))
        jobs.append(thread)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

else:
    PDB_list = []
    AlphaFold_list = []
    for seq in fasta_list:
        PDB(seq, True, config['blast']['PDBdb_path'])
        AmphiProt(seq, True, config['blast']['SwissProtdb_path'])

print(PDB_list)
