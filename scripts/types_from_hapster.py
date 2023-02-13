import sys
import pandas as pd

if __name__ == "__main__":
    plex_map = pd.read_csv(sys.argv[1])
    hapster_dir = sys.argv[2]
    types_outfile = sys.argv[3]

    types_full = pd.DataFrame()
    names = "aliquot,A_1,A_2,B_1,B_2,C_1,C_2,DPA1_1,DPA1_2,DPB1_1,DPB1_2,DQA1_1,DQA1_2,DQB1_1,DQB1_2,DRA_1,DRA_2,DRB1_1,DRB1_2".split(",")
    cases = list(set(plex_map["case"]))
    for case in cases:
        hla_types = pd.read_csv(f"{hapster_dir}/{case}/{case}-DNA-N_haplotype_homozygous.csv", header = None, names = names)
        aliquots = list(set(plex_map[plex_map["case"] == case]["aliquot"]))
        for aliquot in aliquots:
            hla_types["aliquot"] = aliquot
            types_full = pd.concat([types_full, hla_types], ignore_index = True)

    types_full.to_csv(types_outfile, index = False)