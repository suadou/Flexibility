import argparse
from PDB import request
from FASTA import read
import threading


def AmphiProt(query):
    for seq in request.parse_XML(request.BLAST_sp(query)):
        if request.download_AlphaFold(query, seq[0], './files/PDB/') == False:
            continue
        else:
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

print(f"Reading {options.input_file} file.")
fasta_list = read.FASTA_parser(options.input_file)
print(f"{len(fasta_list)} record/s was/were read.")

# Threading for API BLASTp search
threads = len(fasta_list)
jobs = []
for i in range(0, threads):
    thread = threading.Thread(target=AmphiProt(fasta_list[i]))
    jobs.append(thread)

for j in jobs:
    j.start()

for j in jobs:
    j.join()
