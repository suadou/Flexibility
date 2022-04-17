import argparse
from PDB import request
from FASTA import read
from flex_calculus import calculus
import threading
import configparser


# Arguments for the application
parser = argparse.ArgumentParser(
        description="flexibility.py is a standalone application which computes a flexibility score from a protein sequence contained in a FASTA or multi-FASTA file. It is based on experimental data from SwissProt PDB files or models proposed by AlphaFold. It retreats a parse-able text file and a graphical representation of flexibility")
parser.add_argument('-i', '-in', '--input', dest='input_file', action='store',
                    help='FASTA or Multi-FASTA protein sequence input file', required=True, default=None)
parser.add_argument('-o', '-out', '--output', dest='output_file', action='store',
                    help='Protein flexibility output file. It will give 2 files ([].txt - Parseable text file) and ([].png - Graphycal representation of flexibility scores). If 2 or more sequences are provided the names will be []_[FASTA identifier].', required=False, default='flexibility_output')
parser.add_argument('-alpha', '--alpha_threshold', dest='alphafold_threshold', action='store',
                    help='AlphaFold identity threshold. It set the threshold of AlphaFold identity of the query to use it as the true complete structure.', required=False, default=0.95)
options = parser.parse_args()
args, leftovers = parser.parse_known_args()

# Functions defined for the main program

def alphafold(query, local, database=None):
    """
    Calls the function needed to search for match on SwissProt database using BLASTp and download a PDB file of the best match if it exists.
    The function return a PDB object from PDB.request module.
    The arguments are the following:
        query - It is a Query class from FASTA.read module wich contains sequence and identifier.
        local - True for a local search and False for a API-mediated
        database - Path for the local database
    """
    if local == False:
        BLAST = request.blast_sp(query)
        for seq in request.parse_xml(BLAST):
            if request.download_alphafold(query, seq[0], './') == False:
                continue
            else:
                return request.Pdb(
                    query.identifier, seq[0], None, './' + query.identifier + '_AlphaFold.pdb', seq[1], seq[2], seq[3], seq[4])
    else:
        BLAST = request.blast_commandline(query, database)
        BLAST = BLAST[0].split('\n|\r')
        for seq in request.parse_blast(BLAST):
            if request.download_AlphaFold(query, seq[0].split('|')[1], './') == False:
                continue
            else:
                return request.Pdb(
                    query.identifier, seq[0], None, './' + query.identifier + '_AlphaFold.pdb', seq[1], seq[2], seq[3], seq[4])


def pdb(query, local, database=None):
    """
    Calls the function needed to search for match on PDB database using BLASTp and download a PDB file of the best match if it exists.
    The function return a PDB object from PDB.request module.
    The arguments are the following:
        query - It is a Query class from FASTA.read module wich contains sequence and identifier.
        local - True for a local search and False for a API-mediated
        database - Path for the local database
    """
    if local == False:
        BLAST = request.blast_pdb(query)
        for seq in request.parse_xml(BLAST):
            if request.download_pdb(query, seq[0][0:4], './') == False:
                continue
            else:
                return request.Pdb(
                    query.identifier, seq[0], seq[0][5:], './' + query.identifier + '.pdb', seq[1], seq[2], seq[3], seq[4])
    else:
        BLAST = request.blast_commandline(query, database)
        BLAST = BLAST[0].split('\n|\r')
        for seq in request.parse_blast(BLAST):
            if request.download_pdb(query, seq[0][0:4], './') == False:
                continue
            else:
                return request.Pdb(
                    query.identifier, seq[0], seq[0][5:], './' + query.identifier + '.pdb', seq[1], seq[2], seq[3], seq[4])

def retrieving_score(pdb_prefix, chain_id=None):
    """
    Calls the functions needed to calculate the flexibility score of a PDB file. It returns a array with them.
    """
    print(
        f"Computing flexibility score on {pdb_prefix}...")
    matrix = calculus.general_calculation(
        './'+pdb_prefix+'.pdb', chain_id)
    print("Done")
    return matrix



# Read the configure archive to determine local or online BLASTp search
config = configparser.ConfigParser()
config.read('config.ini')
# Read the FASTA file provided
print(f"Reading {options.input_file} file.")
fasta_list = read.fasta_parser(options.input_file)
print(f"{len(fasta_list)} record/s was/were read.")

# Variables to store SwissProt and AlphaFold matches
pdb_list = [[]]
alphafold_list = []


# Threading for API BLASTp. Online version
if config['blast']['local'] == 'False':
    threads = len(fasta_list)
    jobs = []
    # Loop to for all sequences provided
    for i in range(0, threads):
        print(
            f"Searching similar sequences to {fasta_list[i].identifier} in PDB and SwissProt databases using BLASTp server...")
        # Thread to PDB database
        thread = threading.Thread(
            target=pdb_list[i].append(pdb(fasta_list[i], False)))
        jobs.append(thread)
        # Thread to SwissProt database
        thread = threading.Thread(target=alphafold_list.append(
            alphafold(fasta_list[i], False)))
        jobs.append(thread)

    # Start the threading
    for j in jobs:
        j.start()

    # Check if the threading is finished
    for j in jobs:
        j.join()
    
    # For each AlphaFold match check if the identity is higher than -alpha to search for PDB linked to its ids
    i = 0
    for element in alphafold_list:
        if element.identity > float(options.alphafold_threshold):
            pdb_list[i] = []
            print(
                f"The model fits with {alphafold_list[i].identifier} with {alphafold_list[i].identity} of identity \nSearching PDB files linked to the UniProt code")
            request.download_uniprot_xml(
                alphafold_list[i].identifier, './')
            pdb_swissprotid = request.parse_uniprot_xml(
                alphafold_list[i].identifier + '_uniprot.xml')
            if pdb_swissprotid:
                print(f'{len(pdb_swissprotid)} PDB files found')
                for pdb in pdb_swissprotid:
                    print(f"Downloading {pdb[0]} from PDB server...")
                    if request.download_pdb(fasta_list[i], pdb[0], './') == False:
                        print("Cannot download {pdb[0]} PDB file. It was omitted")
                        continue
                    else:
                        pdb_list[i].append(request.Pdb(fasta_list[i].identifier, pdb[0], pdb[1], './'
                                                    + pdb[0] + '_' + fasta_list[i].identifier + '.pdb', 1, pdb[2], pdb[3], []))
                        print("Done")
        i = i + 1

# Local version
elif config['blast']['local'] == 'True':
    # Loop to go for all sequences provided
    for i in range(0, len(fasta_list)):
        # Searching for AlphaFold matches
        print(
            f"Searching similar sequences to {fasta_list[i].identifier} in SwissProt databases using BLASTp...")
        try:
            alphafold_list.append(
                alphafold(fasta_list[i], True, config['blast']['SwissProtdb_path']))
            print(
                f"AlphaFold match: {alphafold_list[i].identifier.split('|')[1]} ({alphafold_list[i].identity:.2f})")
        except IndexError:
            alphafold_list.append(None)
            print("No AlphaFold match was found")
        # Searching for PDB files linked to the SwissProt id if identity is higher than a -alpha threshold
        if alphafold_list[i].identity > float(options.alphafold_threshold):
            print(
                f"The model fits with {alphafold_list[i].identifier.split('|')[1]} with {alphafold_list[i].identity} of identity \nSearching PDB files linked to the UniProt code")
            request.download_uniprot_xml(
                alphafold_list[i].identifier.split('|')[1], './')
            pdb_swissprotid = request.parse_uniprot_xml(
                alphafold_list[i].identifier.split('|')[1] + '_uniprot.xml')
            if pdb_swissprotid:
                print(f'{len(pdb_swissprotid)} PDB files found')
                for pdb in pdb_swissprotid:
                    print(f"Downloading {pdb[0]} from PDB server...")
                    if request.download_pdb(fasta_list[i], pdb[0], './') == False:
                        print("Cannot download {fasta_list[i]} PDB file. It was omitted")
                        continue
                    else:
                        pdb_list[i].append(request.Pdb(fasta_list[i].identifier, pdb[0], pdb[1], './'
                                                    + pdb[0] + '_' + fasta_list[i].identifier + '.pdb', 1, pdb[2], pdb[3], []))
                        print("Done")
            # If no PDB is found linked to Swissprot id, it seeks for PDBs match using BLASTp
            else:
                print(f"No PDB sequences found linked to {alphafold_list[i].identifier.split('|')[1]}")
                try:
                    pdb_list[i].append(
                        pdb(fasta_list[i], True, config['blast']['PDBdb_path']))
                    print(
                        f"PDB match: {pdb_list[i][0].identifier} ({pdb_list[i][0].identity:.2f})")
                except IndexError:
                    pdb_list[i].append(None)
                    print(f"No PDB match was found")    
        # If identity AlphaFold is lower than -alpha it seeks for PDBs matches using BLASTp
        else:
            try:
                pdb_list[i].append(
                    pdb(fasta_list[i], True, config['blast']['PDBdb_path']))
                print(
                    f"PDB match: {pdb_list[i][0].identifier} ({pdb_list[i][0].identity:.2f})")
            except IndexError:
                pdb_list[i].append(None)
                print(f"No PDB match was found")
        if i < len(fasta_list):
            pdb_list.append([])
else:
    raise Exception("Error in config.ini file\nlocal must be \"True\" or \"False\"")

name = options.output_file 
for i in range(0, len(fasta_list)):
    if len(fasta_list) > 1:
	    options.output_file = name + "_" + fasta_list[i].identifier
    plot = True
    matrix_pdb = []
    if alphafold_list[i] != None and pdb_list[i][0] != None:
        if alphafold_list[i].identity > float(options.alphafold_threshold) and len(pdb_list[i]) > 1:
            matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
            matrix_pdb = calculus.general_calculation_multiple(
                pdb_list[i], alphafold_list[i])
        elif alphafold_list[i].identity > pdb_list[i][0].identity:
            matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
        else:
            matrix = retrieving_score(fasta_list[i].identifier,
                             pdb_list[i][0].chain)
    elif alphafold_list[i] != None and pdb_list[i][0] == None:
        matrix = retrieving_score(fasta_list[i].identifier+"_AlphaFold")
    elif alphafold_list[i] == None and pdb_list[i][0] != None:
        matrix = retrieving_score(fasta_list[i].identifier, pdb_list[i][0].chain)
    else:
        print("No PDB files were obtained. Flexibility cannot be calculated")
        plot = False
    if plot:
        print("Saving data and plot...")
        calculus.represent_data(matrix, options.output_file, matrix_pdb)
        print("Done")
