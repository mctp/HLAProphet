# HLAProphet

HLAProphet is a tool that allows for personalized quantification of the HLA proteins in TMT MS/MS data using FragPipe. HLAProphet takes a list of known HLA types for all samples in an experiment and creates a fasta file containing a harmonized list of protein sequences. This HLA fasta file is then appended to an existing reference proteome for use as an augmented search database with FragPipe. After running FragPipe, HLAProphet then uses a modified version of the TMT-integrator algorithm to quantify HLA proteins at the gene and allele level.

# Setup


# HLAProphet workflow part 1:

1) Generate HLA types for all samples using any method, formatted as seen in `examples/types.csv`

2) Create a local version of the IMGT/HLA protein database using `scripts/make_imgt_database.py`
    ```
    python scripts/make_imgt_database.py \
        examples/IMGT \
        examples/IMGT/IMGT_HLA.fa
    ```

3) Create an HLA fasta reference using `scripts/make_hla_fasta.py`. If two separate HLA types produce the same protein product, the protein is only included once in the output database. A relationship table is produced to tie original HLA types to the matching sequence in the HLA fasta, after clashes are resolved.
    ```
    python scripts/make_hla_fasta.py \
        examples/types.csv \
        examples/IMGT/IMGT_HLA.fa \
        examples/example_HLA.fa \
        examples/example_relationships.csv
    ```

3) Predict tryptic peptides using `scripts/tryptic_peptides.py`
    ```
    python scripts/tryptic_peptides.py \
        examples/example_HLA.fa \
        2 \
        7 \
        50 \
        examples/example_tryptic_peptides.csv
    ```

# Fragpipe workflow
    1. Create personalized database using philosopher. This step combines a standard reference proteome (i.e. GENCODE)
    with the cohort personalized HLA reference produced by HLAProphet.

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