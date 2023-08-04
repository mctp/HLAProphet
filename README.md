# HLAProphet

HLAProphet is a tool that allows for personalized quantification of the HLA proteins in TMT MS/MS data using FragPipe. HLAProphet takes a list of known HLA types for all samples in an experiment and creates a fasta file containing a harmonized list of protein sequences. This HLA fasta file is then appended to an existing reference proteome for use as an augmented search database with FragPipe. After running FragPipe, HLAProphet then uses a modified version of the TMT-integrator algorithm to quantify HLA proteins at the gene and allele level.

# Setup

1) Install the HLAProphet conda environment
    ```
    conda install --file HLAProphet.yml
    ```
2) Install FragPipe: https://fragpipe.nesvilab.org/

3) Install MSFragger: https://msfragger.nesvilab.org/

4) Install Philosopher: https://philosopher.nesvilab.org/

# HLAProphet workflow part 1:

1) Generate HLA types for all samples using any method, formatted as seen in `examples/types.csv`

2) Create a local version of the IMGT/HLA protein database using `scripts/make_imgt_database.py`
    ```
    #DOWNLOAD_DIR is a folder where the source IMGT reference files will be downloaded
    #IMGT_FASTA is the filename of the final combined fasta file containing all HLA sequences
    python scripts/make_imgt_database.py \
        $DOWNLOAD_DIR \
        $IMGT_FASTA
    ```

3) Create an HLA fasta reference using `scripts/make_hla_fasta.py`. If two separate HLA types produce the same protein product, the protein is only included once in the output database. A relationship table is produced to tie original HLA types to the matching sequence in the HLA fasta, after clashes are resolved.
    ```
    #HLA_TYPES is a table containing HLA types of all samples in an experiment
    #IMGT_FASTA is the combined IMGT database fasta file created in the previous step
    #HLA_FASTA is the output filename for the fasta file containing all HLA sequences for the experiment
    #HLA_RELATIONSHIPS is a relationship table matching original HLA types in $HLA_TYPES to the condensed HLA types in $HLA_FASTA
    python scripts/make_hla_fasta.py \
        $HLA_TYPES \
        $IMGT_FASTA \
        $HLA_FASTA \
        $HLA_RELATIONSHIPS
    ```

3) Predict tryptic peptides using `scripts/tryptic_peptides.py`
    ```
    #HLA_FASTA is the HLA fasta file created in the previous step
    #N_MISSED_CLEAVAGE is the number of allowed missed cleavages by trypsin, suggested value is 2
    #MIN_PEPTIDE_LENGTH is the minimum tryptic peptide length to keep, suggested value is 7
    #MAX_PEPTIDE_LENGTH is the maximum tryptic peptide length to keep, suggested value is 50
    #TRYPTIC_PEPTIDES is the filename of the final output file containing all predicted tryptic peptides for all HLA sequences
    python scripts/tryptic_peptides.py \
        $HLA_FASTA \
        $N_MISSED_CLEAVAGE \
        $MIN_PEPTIDE_LENGTH \
        $MAX_PEPTIDE_LENGTH \
        $TRYPTIC_PEPTIDES
    ```

# Fragpipe workflow
1. Create personalized database using philosopher. This step combines a standard reference proteome (i.e. GENCODE) with the cohort personalized HLA reference produced by HLAProphet. Be sure to remove existing HLA sequences in the reference proteome before combining.

    ```
    cd examples
    philosopher workspace --init
    philosopher database --custom $GENCODE --add example_HLA.fa --contam
    philosopher workspace --clean
    cd ../
    ```

2. Save sample manifest in FragPipe using the GUI.

3. Save experiment workflow in FragPipe using the GUI. 

4. Run FragPipe in headless mode

    ```
    fragpipe --headless \
        --manifest $MANIFEST \
        --workflow $WORKFLOW \
        --threads $NCORES \
        --workdir $FRAGPIPE_WORKDIR \
        --config-msfragger $MSFRAGGER \
        --config-philosopher $PHILOSOPHER
    ```

# HLAProphet workflow part 2:
1) Run HLA quant.

    ```
    python scripts/hla_quant.py \
        examples/example_relationships.csv \
        examples/example_tryptic_peptides.csv \
        $FRAGPIPE_WORKDIR \
        $REF_PATTERN \
        $PLEX_SIZE \
        $POOL_N \
        $OUTDIR  \
        $OUT_PREFIX
    ```
