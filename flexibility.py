import argparse
from PDB import request
from FASTA import read
from flex_calculus import calculus
import threading
import configparser


def AlphaFold(query, local, database=None):
    if local == False:
        BLAST = request.BLAST_sp(query)
        for seq in request.parse_XML(BLAST):
            if request.download_AlphaFold(query, seq[0].split('|')[1], './files/PDB/') == False:
                continue
            else:
                return request.PDB(
                    query.identifier, seq[0], None, './files/PDB/' + query.identifier + '_AlphaFold.pdb', seq[1], seq[2], seq[3], seq[4])
    else:
        BLAST = request.BLAST_commandline(query, database)
        BLAST = BLAST[0].split('\n|\r')
        for seq in request.parse_blast(BLAST):
            if request.download_AlphaFold(query, seq[0].split('|')[1], './files/PDB/') == False:
                continue
            else:
                return request.PDB(
                    query.identifier, seq[0], None, './files/PDB/' + query.identifier + '_AlphaFold.pdb', seq[1], seq[2], seq[3], seq[4])


def PDB(query, local, database=None):
    if local == False:
        BLAST = request.BLAST_PDB(query)
        for seq in request.parse_XML(BLAST):
            if request.download_PDB(query, seq[0][0:4], './files/PDB/') == False:
                continue
            else:
                return request.PDB(
                    query.identifier, seq[0], seq[0][5:], './files/PDB/' + query.identifier + '.pdb', seq[1], seq[2], seq[3], seq[4])
    else:
        BLAST = request.BLAST_commandline(query, database)
        BLAST = BLAST[0].split('\n|\r')
        for seq in request.parse_blast(BLAST):
            if request.download_PDB(query, seq[0][0:4], './files/PDB/') == False:
                continue
            else:
                return request.PDB(
                    query.identifier, seq[0], seq[0][5:], './files/PDB/' + query.identifier + '.pdb', seq[1], seq[2], seq[3], seq[4])


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
args, leftovers = parser.parse_known_args()

config = configparser.ConfigParser()
config.read('config.ini')
print(f"Reading {options.input_file} file.")
fasta_list = read.FASTA_parser(options.input_file)
print(f"{len(fasta_list)} record/s was/were read.")
PDB_list = [[]]
AlphaFold_list = []
# Threading for API BLASTp search
if config['blast']['local'] == 'False':
    threads = len(fasta_list)
    jobs = []
    for i in range(0, threads):
        print(
            f"Searching similar sequences to {fasta_list[i].identifier} in PDB and SwissProt databases using BLASTp server...")
        thread = threading.Thread(
            target=PDB_list.append(PDB(fasta_list[i], False)))
        jobs.append(thread)
        thread = threading.Thread(target=AlphaFold_list.append(
            AlphaFold(fasta_list[i], False)))
        jobs.append(thread)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

else:
    for i in range(0, len(fasta_list)):
        print(
            f"Searching similar sequences to {fasta_list[i].identifier} in SwissProt databases using BLASTp...")
        try:
            AlphaFold_list.append(
                AlphaFold(fasta_list[i], True, config['blast']['SwissProtdb_path']))
            print(
                f"AlphaFold match: {AlphaFold_list[i].identifier.split('|')[1]} ({AlphaFold_list[i].identity:.2f})")
        except IndexError:
            AlphaFold_list.append(None)
            print("No AlphaFold match was found")
        if AlphaFold_list[i].identity > 0.95:
            print(
                f"The model fits with {AlphaFold_list[i].identifier.split('|')[1]} with {AlphaFold_list[i].identity} of identity \nSearching PDB files linked to the UniProt code")
            request.download_uniprot_xml(
                AlphaFold_list[i].identifier.split('|')[1], './')
            pdb_list = request.parse_uniprot_xml(
                AlphaFold_list[i].identifier.split('|')[1] + '_uniprot.xml')
            if pdb_list:
                print(f'{len(pdb_list)} PDB files found')
                for pdb in pdb_list:
                    print(f"Downloading {pdb[0]} from PDB server...")
                    if request.download_PDB(fasta_list[i], pdb[0], './') == False:
                        print("Cannot download {fasta_list[i]} PDB file. It was omitted")
                        continue
                    else:
                        PDB_list[i].append(request.PDB(fasta_list[i].identifier, pdb[0], pdb[1], './'
                                                    + pdb[0] + '_' + fasta_list[i].identifier + '.pdb', 1, pdb[2], pdb[3], []))
                        print("Done")
            else:
                try:
                    PDB_list.append(
                        PDB(fasta_list[i], True, config['blast']['PDBdb_path']))
                    print(
                        f"PDB match: {PDB_list[i].identifier} ({PDB_list[i].identity:.2f})")
                except IndexError:
                    PDB_list.append(None)
                    print(f"No PDB match was found")    
        else:
            try:
                PDB_list.append(
                    PDB(fasta_list[i], True, config['blast']['PDBdb_path']))
                print(
                    f"PDB match: {PDB_list[i].identifier} ({PDB_list[i].identity:.2f})")
            except IndexError:
                PDB_list.append(None)
                print(f"No PDB match was found")


def retrieving_score(pdb_prefix, chain_id=None):
    print(
        f"Computing flexibility score on {pdb_prefix}...")
    matrix = calculus.general_calculation(
        './files/PDB/'+pdb_prefix+'.pdb', chain_id)
    print("Done")
    return matrix


for i in range(0, len(fasta_list)):
    plot = True
    if AlphaFold_list[i] != None and PDB_list[i] != None:
        if AlphaFold_list[i].identity > 0.95 and pdb_list:
            matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
            matrix_pdb = calculus.general_calculation_multiple(
                PDB_list[i], AlphaFold_list[i])
        elif AlphaFold_list[i].identity > PDB_list[i].identity:
            matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
        else:
            matrix = retrieving_score(fasta_list[i].identifier,
                             PDB_list[i].chain)
    elif AlphaFold_list[i] != None and PDB_list[i] == None:
        matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
    elif AlphaFold_list[i] == None and PDB_list[i] != None:
        matrix = retrieving_score(fasta_list[i].identifier, PDB_list[i].chain)
    else:
        print("No PDB files were obtained. Flexibility cannot be calculated")
        plot = False
    if plot:
        print("Saving data and plot...")
        calculus.represent_data(matrix, "./", matrix_pdb)
        print("Done")
