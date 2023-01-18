# HLAProphet workflow part 1:
## TODO: Make these steps a single command from HLAProphet that takes as input HLA types

    1) Generate HLA types for all samples using any method (Hapster used for this example), 
    formatted as seen in LSCC_types.csv

    2) Create an HLA fasta reference using make_hla_fasta.py

    3) Predict tryptic peptides using tryptic_peptides.py

# Fragpipe workflow
    1. Create personalized database using philosopher. This step combines gencode (HLA seqs removed)
    with the cohort personalized HLA reference produced by HLAProphet.

    ```
    philosopher workspace --init
    philosopher database --custom $GENCODE --add $HLA --contam
    philosopher workspace --clean
    ```

    2. Save sample manifest in FragPipe using the GUI.

    3. Save experiment workflow in FragPipe using the GUI. 

    4. Run FragPipe in headless mode

    ```
    COHORT=LSCC
    MANIFEST=/mctp/share/users/mumphrey/HLAProphet/data/experiments/$COHORT/ms/Phosphoproteome/${COHORT}.HLAProphet.PhosphoGlyco.manifest
    WORKFLOW=/mctp/share/users/mumphrey/HLAProphet/data/experiments/$COHORT/ms/Phosphoproteome/${COHORT}.HLAProphet.PhosphoGlyco.workflow
    WORKDIR=/mctp/share/users/mumphrey/HLAProphet/data/experiments/${COHORT}/ms/Phosphoproteome
    MSFRAGGER=/mctp/share/users/mumphrey/HLAProphet/bin/MSFragger-3.6/MSFragger-3.6.jar
    PHILOSOPHER=/mctp/share/users/mumphrey/HLAProphet/bin/philosopher_4.6
    fragpipe --headless \
        --manifest $MANIFEST \
        --workflow $WORKFLOW \
        --threads 32 \
        --workdir $WORKDIR \
        --config-msfragger $MSFRAGGER \
        --config-philosopher $PHILOSOPHER \
        --dry-run
    ```

# HLAProphet workflow part 2:

    1) 