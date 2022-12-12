import pandas as pd
import numpy as np
import glob
import sys

def hla_from_psm(filename, n):
    print(f"Processing: {filename}")
    #Apply initial filters across all PSMs
    print(f"\tReading PSM")
    psm = pd.read_csv(filename, sep = "\t")
    #Only keep PSMs with Ref intensity > 0, Purity > 0.5, PeptidePropeht Prob  > 0.5
    print(f"\tFiltering PSM")
    psm = psm[(psm[ref_name] != 0) & (psm["Purity"] > 0.5) & (psm["PeptideProphet Probability"] > 0.5)]
    #Drop the bottom 5th percentile of MS2 intensities
    psm["MS2"] = psm.iloc[:,-plex_size:].sum(1)
    psm = psm[psm["MS2"] > np.quantile(psm["MS2"], q = 0.05)]
    #Keep only the highest MS2 PSM for each Peptide
    psm = psm.loc[psm.groupby("Peptide")["MS2"].idxmax()]
    #After filtering, move to just HLA PSMs
    #Discard any PSMs that have potential matches outside the HLAs of interest
    cols = ["Peptide", "Protein Start", "Protein End", "Protein", "Mapped Proteins", "Intensity"] + psm.iloc[:, -plex_size - 1:-1].columns.tolist()
    psm_hla = psm[psm["Protein"].str.startswith("HLA")].loc[:, cols].copy()
    psm_hla = pd.melt(psm_hla, id_vars = cols[:6], value_vars = cols[6:], var_name = "Aliquot", value_name = "MS2")
    psm_hla["Mapped Proteins"] = psm_hla["Mapped Proteins"].fillna("")
    psm_hla["Proteins"] = psm_hla[["Protein", "Mapped Proteins"]].apply(", ".join, axis = 1)
    psm_hla = psm_hla[["Peptide", "Protein Start", "Protein End", "Aliquot", "Proteins", "Intensity", "MS2"]]
    #For each PSM in each aliquot, check if the peptide is predicted to be present in the sample
    psm_hla["Predicted"] = False
    psm_hla["Predicted_n"] = 0
    psm_hla["Aliquot_prot"] = [[] for _ in range(len(psm_hla))]
    print(f"\tChecking peptide predictions")
    for i, row in psm_hla.iterrows():
        aliquot_alleles = hla_types[hla_types["aliquot"] == row["Aliquot"]]["adjusted"]
        for allele in aliquot_alleles:
            if row["Proteins"].find(allele) != -1:
                psm_hla.loc[i, "Predicted"] = True
                psm_hla.loc[i, "Predicted_n"] += 1
                psm_hla.loc[i, "Aliquot_prot"].append(allele)
    #Record plex number for later use
    psm_hla["Plex"] = n
    return psm_hla        

def main():
    fragpipe_workdir = sys.argv[1]
    ref_name = sys.argv[2]
    plex_size = sys.argv[3]
    hla_types = pd.read_csv(sys.argv[4])
    tryptic_predictions = pd.read_csv(sys.argv[5])
    outfile_prefix = sys.argv[6]

    hla_peptides = pd.DataFrame()
    psm_filenames = sorted(glob.glob(fragpipe_workdir + "/*/psm.tsv"))
    for i, f in enumerate(psm_filenames):
        hla_peptides = pd.concat([hla_peptides, hla_from_psm(f, i)])

if __name__ == "__main__":
    main()
