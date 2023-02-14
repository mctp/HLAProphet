import sys
import subprocess
from Bio import SeqIO

if __name__ == "__main__":
    tmp_dir = sys.argv[1]
    output_filename = sys.argv[2]

    #Make sure output locations exist
    pathlib.Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_filename).parent.mkdir(parents=True, exist_ok=True)

    genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1"]

    for gene in genes:
        try:
            subprocess.run(f"wget -O {tmp_dir}/{gene}.fa https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/{gene}_prot.fasta".split(" "))
        except:
            f"Can't download fasta for HLA-{gene}"

    imgt_database = []
    seen = []
    for gene in genes:
        tmp_fasta = [x for x in SeqIO.parse(f"{tmp_dir}/{gene}.fa", 'fasta')]
        for record in tmp_fasta:
            type_full = record.description.split(" ")[1]
            type_short = ":".join(type_full.split(":")[:2])
            if not type_short in seen:
                seen.append(type_short)
                record.id = type_short
                imgt_database.append(record)
    
    with open(output_filename, 'w') as of:
        SeqIO.write(imgt_database, of, "fasta")
        