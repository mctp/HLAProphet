# HLAProphet workflow part 1:
## TODO: Make these steps a single command from HLAProphet that takes as input HLA types

1) Generate HLA types for all samples using any method (Hapster used for this example), 
formatted as seen in LSCC_types.csv
    ```
    python scripts/types_from_hapster.py \
        data/experiments/$COHORT/raw/${COHORT}_plex_map.csv \
        $HAPSTER_DIR \
        data/experiments/$COHORT/raw/${COHORT}_types.csv
    ```

2) Create an HLA fasta reference using make_hla_fasta.py
    ```
    python scripts/make_hla_fasta.py \
        data/experiments/$COHORT/raw/${COHORT}_types.csv \
        data/IMGT/IMGT_HLA.fa \
        data/experiments/$COHORT/refs/${COHORT}_HLA.fa \
        data/experiments/$COHORT/refs/${COHORT}_relationships.csv
    ```

3) Predict tryptic peptides using tryptic_peptides.py
    ```
    python scripts/tryptic_peptides.py \
        data/experiments/$COHORT/refs/${COHORT}_HLA.fa \
        2 \
        7 \
        50 \
        data/experiments/$COHORT/refs/${COHORT}_tryptic_peptides.csv
    ```

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
    EXPERIMENT=Proteome
    MANIFEST=data/experiments/$COHORT/ms/$EXPERIMENT/${COHORT}.HLAProphet.manifest
    WORKFLOW=data/experiments/$COHORT/ms/$EXPERIMENT/${COHORT}.HLAProphet.workflow
    WORKDIR=data/experiments/${COHORT}/ms/$EXPERIMENT
    MSFRAGGER=bin/MSFragger-3.7/MSFragger-3.7.jar
    PHILOSOPHER=bin/philosopher
    NCORES=32
    fragpipe --headless \
        --manifest $MANIFEST \
        --workflow $WORKFLOW \
        --threads $NCORES \
        --workdir $WORKDIR \
        --config-msfragger $MSFRAGGER \
        --config-philosopher $PHILOSOPHER
    ```

# HLAProphet workflow part 2:
    1) Run HLA quant

    ```
    python scripts/hla_quant.py \
        data/experiments/$COHORT/refs/${COHORT}_relationships.csv \
        data/experiments/$COHORT/refs/${COHORT}_tryptic_peptides.csv \
        data/experiments/$COHORT/ms/$EXPERIMENT \
        $REF_PATTERN \
        $PLEX_SIZE \
        $POOL_N \
        data/experiments/$COHORT/tables/$EXPERIMENT \
        $COHORT
    ```