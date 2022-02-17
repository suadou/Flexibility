import Bio
import argparse

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
option = parser.parse_args()
