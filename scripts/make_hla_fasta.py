"""
Script takes a list of sample level HLA types and produces a fasta database
with all observed protein sequences.

Prevents duplicate records in fasta database and creates a relationship
table for protein -> case inside the plex.
"""

import sys
import pandas as pd
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

if __name__ == "__main__":
    hla_types = pd.read_csv(sys.argv[1])
    imgt_database = {record.id:record for record in SeqIO.parse(sys.argv[2], "fasta")}
    fasta_database_filename = sys.argv[3]
    relationship_database_filename = sys.argv[4]

    fasta_database = {}
    relationship_database = {aliq:{"original":[], "adjusted":[]} for aliq in hla_types["aliquot"]}

    genes = hla_types.columns[2:]
    for aliq in hla_types["aliquot"]:
        tmp_types = [hla_types[hla_types["aliquot"] == aliq][g].values[0] for g in genes]
        for t in tmp_types:
            if isinstance(t, float):
                continue
            t_short = ":".join(t.split(":")[:2])
            seq = imgt_database[t_short]
            if seq.seq in fasta_database.values():
                #If HLA sequence is present under a different name, just let all
                #samples use the first observed name.
                existing_name = [k for k, v in fasta_database.items() if v == seq.seq][0]
                relationship_database[aliq]["adjusted"].append(existing_name)
                relationship_database[aliq]["original"].append(t)
            else:
                #Else it is a new sequence
                relationship_database[aliq]["adjusted"].append(seq.id)
                relationship_database[aliq]["original"].append(t)
                fasta_database[seq.id] = seq.seq

    with open(fasta_database_filename, 'w') as of:
        fasta_database_records = []
        i = 0
        for seq_id, seq in sorted(fasta_database.items()):
            fasta_database_records.append(SeqRecord(id = f"HLA-{seq_id}".replace("*", "-").replace(":", "-"), seq = seq, description = ""))
            i += 1
        SeqIO.write(fasta_database_records, of, "fasta")

    with open(relationship_database_filename, 'w') as of:
        of.write(f"aliquot,gene,original,adjusted\n")
        for case in relationship_database.keys():
            original_alleles = relationship_database[case]["original"]
            adjusted_alleles = relationship_database[case]["adjusted"]
            for original, adjusted in zip(original_alleles, adjusted_alleles):
                gene = adjusted.split("*")[0]
                of.write(f"{case},HLA-{gene},HLA-{original},HLA-{adjusted}\n")
