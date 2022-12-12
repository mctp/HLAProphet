import sys
from Bio import SeqIO
from pyteomics import parser

if __name__ == "__main__":
    fasta = [x for x in SeqIO.parse(sys.argv[1], "fasta")]
    missed = int(sys.argv[2])
    min_len = int(sys.argv[3])
    max_len = int(sys.argv[4])
    outfile = sys.argv[5]

    with open(outfile, 'w') as of:
        of.write(f"allele,peptide\n")
        for record in fasta:
            allele = record.id
            peptides = parser.cleave(sequence = str(record.seq), rule = parser.expasy_rules["trypsin"], missed_cleavages = missed, min_length = min_len, max_length = max_len)
            for peptide in peptides:
                of.write(f"{allele},{peptide}\n")