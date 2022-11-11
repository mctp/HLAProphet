"""
Script takes a list of sample level HLA types and produces a fasta database
with all observed protein sequences.

Prevents duplicate records in fasta database and creates a relationship
table for protein -> case inside the plex.
"""

import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

if __name__ == "__main__":
    hla_types = pd.read_csv(sys.argv[1])
    imgt_database = {record.id:record for record in SeqIO.parse(sys.argv[2], "fasta")}
    plex_database_filename = sys.argv[3]
    relationship_database_filename = sys.argv[4]

    plex_database = {}
    relationship_database = {aliq:{"original":[], "adjusted":[]} for aliq in hla_types["aliquot"]}

    genes = hla_types.columns[2:]
    for aliq in hla_types["aliquot"]:
        tmp_types = [hla_types[hla_types["aliquot"] == aliq][g].values[0] for g in genes]
        for t in tmp_types:
            t_short = ":".join(t.split(":")[:2])
            seq = imgt_database[t_short]
            if seq.seq in plex_database.values():
                #If HLA sequence is present under a different name, just let all
                #samples use the first observed name.
                existing_name = [k for k, v in plex_database.items() if v == seq.seq][0]
                relationship_database[aliq]["adjusted"].append(existing_name)
                relationship_database[aliq]["original"].append(t)
            else:
                #Else it is a new sequence
                relationship_database[aliq]["adjusted"].append(seq.id)
                relationship_database[aliq]["original"].append(t)
                plex_database[seq.id] = seq.seq

    with open(plex_database_filename, 'w') as of:
        plex_database_records = []
        i = 0
        for seq_id, seq in sorted(plex_database.items()):
            #Add in dummy fields corresponding to what would normally be uniprot IDs:
            #>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
            #Replace * and : characters to prevent problems parsing nonstandard characters
            #Unique ID has the format [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
            unique_id = f"Q0HLA0Z{i:03}"
            #dummy_id = f"sp|{unique_id}|HLA-{seq_id} HLA-{seq_id} OS=Homo Sapiens GN=HLA-{seq_id} PE=2 SV=1".replace("*", "-").replace(":", "-")
            #dummy_id = f"{i+1}|{i+2}|{i+3}|{i+4}|{i+5}|HLA-{seq_id}|HLA-{seq_id}|{len(seq)}".replace("*", "-").replace(":", "-")
            dummy_id = f"HLA-{seq_id}".replace("*", "-").replace(":", "-")
            plex_database_records.append(SeqRecord(id = dummy_id, seq = seq, description = ""))
            i += 1
        SeqIO.write(plex_database_records, of, "fasta")

    with open(relationship_database_filename, 'w') as of:
        of.write(f"case,gene,original,adjusted\n")
        for case in relationship_database.keys():
            original_alleles = relationship_database[case]["original"]
            adjusted_alleles = relationship_database[case]["adjusted"]
            for original, adjusted in zip(original_alleles, adjusted_alleles):
                gene = adjusted.split("*")[0]
                of.write(f"{case},HLA-{gene},HLA-{original},HLA-{adjusted}\n")
