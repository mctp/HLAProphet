import pandas as pd
import numpy as np
import glob
import sys
import scipy

def parse_psm(filename, plex_size, ref_name, hla_types, tryptic_predictions, n):
    #Parses PSMs to get HLA peptide intensities, and global aliquot ratios for later quantification
    print(f"Processing: {filename}")
    #Apply initial filters across all PSMs
    print(f"\tReading PSM")
    psm = pd.read_csv(filename, sep = "\t")
    samples = [x for x in psm.columns.values[-plex_size:] if x != ref_name]
    #Only keep PSMs with Ref intensity > 0, Purity > 0.5, PeptidePropeht Prob  > 0.5
    print(f"\tFiltering PSM")
    psm = psm[(psm[ref_name] != 0) & (psm["Purity"] > 0.5) & (psm["PeptideProphet Probability"] > 0.5)]
    #Drop the bottom 5th percentile of MS2 intensities
    psm["MS2"] = psm.iloc[:,-plex_size:].sum(1)
    psm = psm[psm["MS2"] > np.quantile(psm["MS2"], q = 0.05)]
    #Keep only the highest MS2 PSM for each Peptide
    psm = psm.loc[psm.groupby("Peptide")["MS2"].idxmax()]
    #Record ratios for all samples
    ratios = psm[samples].div(psm[ref_name], axis = 0).median(axis = 0).to_frame(name = "med_ratio").reset_index(names = ["Aliquot"])
    #Parse HLA PSMs
    cols = ["Peptide", "Protein Start", "Protein End", "Protein", "Mapped Proteins", "Intensity"] + samples + [ref_name]
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
        predicted_types = set(tryptic_predictions[tryptic_predictions["peptide"] == row["Peptide"]]["allele"])
        for allele in aliquot_alleles:
            if allele in predicted_types:
                psm_hla.loc[i, "Predicted"] = True
                psm_hla.loc[i, "Predicted_n"] += 1
                psm_hla.loc[i, "Aliquot_prot"].append(allele)
    #Record plex number for later use
    psm_hla["Plex"] = n
    return psm_hla, ratios

def remove_bad_peptides(peptides):
    print(f"Removing peptides with bad signal to noise ratio")
    #Calculate noise as the mean MS2 intensity of all peptides that are not predicted to be present in a sample
    #Don't include in calculations aliquots for which no typing is known
    predicted_aliquots = set(peptides[peptides["Predicted"] == True]["Aliquot"])
    noise_peptides = peptides[(peptides["Predicted"] == False) & (peptides["Aliquot"].isin(predicted_aliquots))]
    noise = noise_peptides.groupby(["Peptide", "Plex"])["MS2"].apply(np.mean).reset_index().rename(columns = {"MS2":"Noise"})
    peptides = peptides.merge(noise, how = "left")
    peptides["Noise_LR"] = np.log2(peptides["MS2"]/peptides["Noise"])
    #Get the sample of all logratios for peptides that are predicted to be absent
    #Ignore cases with no typing information, and cases with MS2 of 0
    noise_ratios = peptides[(peptides["Aliquot"].isin(predicted_aliquots)) & (peptides["Predicted"] == False) & (~peptides["Noise_LR"].isna()) & (np.isfinite(peptides["Noise_LR"]))]["Noise_LR"]
    #For each peptide with an observed noise value, see if predicted peptide LRs
    #appear to come from a different distribution than noise LRs
    ks_tests = pd.DataFrame()
    for peptide in noise["Peptide"].unique():
        signal_ratios = peptides[(peptides["Peptide"] == peptide) & (peptides["Predicted"] == True) & (~peptides["Noise_LR"].isna())]["Noise_LR"]
        if len(signal_ratios) == 0:
            continue
        res = scipy.stats.kstest(signal_ratios, noise_ratios)[1]
        ks_tests = pd.concat([ks_tests, pd.DataFrame({"Peptide":[peptide], "ks_p":[res]})], ignore_index = True)
    return peptides[~peptides["Peptide"].isin(ks_tests[ks_tests["ks_p"] > .05]["Peptide"])]

def main():
    fragpipe_workdir = sys.argv[1]
    ref_name = sys.argv[2]
    plex_size = int(sys.argv[3])
    hla_types = pd.read_csv(sys.argv[4])
    tryptic_predictions = pd.read_csv(sys.argv[5])
    outfile_prefix = sys.argv[6]

    #Get HLA peptide intensities from PSMs, and keep track of global ratio medians for later normalization
    hla_peptides = pd.DataFrame()
    aliquot_ratios = pd.DataFrame()
    psm_filenames = sorted(glob.glob(fragpipe_workdir + "/*/psm.tsv"))
    for i, f in enumerate(psm_filenames):
        peptides, ratios = parse_psm(f, plex_size, ref_name, hla_types, tryptic_predictions, i + 1)
        hla_peptides = pd.concat([hla_peptides, peptides], ignore_index = True)
        aliquot_ratios = pd.concat([aliquot_ratios, ratios], ignore_index = True)
    print(hla_peptides)
    print(aliquot_ratios)

    #Remove peptides that don't have a strong signal to noise ratio
    hla_peptides = remove_bad_peptides(hla_peptides)

    

if __name__ == "__main__":
    main()
