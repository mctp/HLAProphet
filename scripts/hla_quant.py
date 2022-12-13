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
    ratios["Plex"] = n
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
    #Only keep peptides that pass the ks test
    peptides_keep = peptides[~peptides["Peptide"].isin(ks_tests[ks_tests["ks_p"] > .05]["Peptide"])]
    #Only keep peptides that are predicted to be present
    peptides_keep = peptides_keep[peptides_keep["Predicted"]]
    return peptides_keep.drop(columns = ["Noise", "Noise_LR"])

def adjusted_ratio(ratios, aliquot_ratios, pool_n):
    #Unlike most peptides, HLA peptides are present at variable numbers of copies
    #in the reference channel. Peptides that are rare in the population may only be
    #in a couple of samples, so when pooling in the reference channel these peptides
    #become highly diluted. Dilution causes reference channel signal to be low, which
    #causes the intensity ratio to be artificially high due to a smaller denominator.
    #We adjust each ratio taking into account the expected copy number in the reference
    #to reverse the effect of dilution.
    ratios["Ratio"] = ratios["MS2"] / ratios["RefMS2"]
    ratios = ratios.merge(aliquot_ratios, how = "left")
    #First count how many copies we expect to be in the reference channel
    #This result will be imperfect if the HLA type of all samples is not known
    pep_n = ratios[["Peptide", "Aliquot", "Predicted_n"]].drop_duplicates().groupby("Peptide").sum("Predicted_n").reset_index().rename(columns = {"Predicted_n":"Total_n"})
    ratios = ratios.merge(pep_n, how = "left")
    #The direct adjustment would be Ratio * (Peptide_n / Aliquot_n) to reverse the effect of dilution,
    #however this gets too extreme as Peptide_n approaches 0. We therefore add a constant C to both the numerator
    #and denominator to stabilize the adjustment. We will try all values of C 1:10000 to get the value that most reduces
    #the relationship between copy number and ratio
    c_best = 0
    for c_cur in range(10000):
        if c_cur % 100 == 0:
            print(f"Testing C value: {c_cur}")
        c_cur = 0
        ratio_adj = ratios["Ratio"] * ((ratios["Total_n"] + c_cur)/(pool_n + c_cur))
        #Perform our tests on median normalized ratios
        ratio_adj_md = np.log2(ratio_adj) - 
        res = scipy.stats.linregress(ratio_adj, ratios["Total_n"])


    return peptides_adjusted


def calc_ratios(peptides, ref_name, aliquot_ratios, pool_n):
    #For each peptide, take the ratio of the MS2 intensity relative to the reference channel
    refs = peptides[peptides["Aliquot"] == ref_name][["Peptide", "Plex", "MS2"]].copy().rename(columns = {"MS2":"RefMS2"})
    ratios = peptides[peptides["Aliquot"] != ref_name].copy().merge(refs, how = "left")
    ratios_adjusted = adjusted_ratio(ratios, aliquot_ratios, pool_n)

    peptides["LR"] = np.log2(peptides["Ratio"])
    #Also provide median normalized ratios
    
    peptides["LR_MD"] = peptides["LR"] - np.log2(peptides["med_ratio"])
    peptides["Ratio_MD"] = 2**peptides["LR_MD"]
    return peptides[["Peptide", "Protein Start", "Protein End", "Aliquot", "Proteins", "Intensity",
                     "MS2", "Predicted", "Predicted_n", "Aliquot_prot", "Plex", "RefMS2", "Ratio",
                     "Ratio_MD", "LR", "LR_MD"]]

def main():
    hla_types = pd.read_csv(sys.argv[1]) #Relationship database constructed by make_hla_fasta.py
    tryptic_predictions = pd.read_csv(sys.argv[2]) #Tryptic predictions made by tryptic_peptides.py
    fragpipe_workdir = sys.argv[3] #Directory that Fragpipe was fun on
    ref_name = sys.argv[4] #Name of the ref channel in each psm.tsv
    plex_size = int(sys.argv[5]) #Number of channels in each plex for this experiment
    pool_n = int(sys.argv[6]) #Number of aliquots constrbuting to the ref pool for this experiment
    outdir = sys.argv[7]
    outfile_prefix = sys.argv[8] #Prefix for all output files

    #Get HLA peptide intensities from PSMs, and keep track of global ratio medians for later normalization
    hla_peptides = pd.DataFrame()
    aliquot_ratios = pd.DataFrame()
    psm_filenames = sorted(glob.glob(fragpipe_workdir + "/*/psm.tsv"))
    for i, f in enumerate(psm_filenames):
        peptides, ratios = parse_psm(f, plex_size, ref_name, hla_types, tryptic_predictions, i + 1)
        hla_peptides = pd.concat([hla_peptides, peptides], ignore_index = True)
        aliquot_ratios = pd.concat([aliquot_ratios, ratios], ignore_index = True)

    #Remove peptides that don't have a strong signal to noise ratio
    hla_peptides_filtered = remove_bad_peptides(hla_peptides)

    #Calculate peptide ratios relative to the ref channel
    hla_ratios = calc_ratios(hla_peptides_filtered, ref_name, aliquot_ratios, pool_n)
    hla_ratios.to_csv(f"{outdir}/{outfile_prefix}_peptide_ratios.csv", index = False)

if __name__ == "__main__":
    main()
