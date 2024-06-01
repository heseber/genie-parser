#!/usr/bin/env python

##### THIS SCRIPT NEEDS >32GB <64GB of RAM AND RUNS FOR ABOUT 2 MINUTES #####

import os

import genie.genie as gd
import pandas as pd

#########################################################################
# Get the mutation classifications from colleague's file
#########################################################################

# This file defines different ways how to filter mutations when defining the
# gene level status of HER2. Here we read this table.
#
# curated:   TRUE if variant has been curated as activating HER2 mutation
# phase3:    TRUE if variant has been selected as eligable for inclusion in NSCLC
#            Phase 3 clinical study
# tkd:       TRUE if variant is located in the Tyrosine Kinase Domain of HER2;
#            TKD is defined as exons 18 - 22; Pfam TKD modle matches AAs 721 - 975
# exon20ins: TRUE, if variant is an inframe exon20 insertion mutation
# missense:  TRUE, if variant is a SNV missense mutation

mut_classes = (
    pd.read_table("from_her2_activating_mutations.tsv", comment="#")
    .rename(columns={"aa_change": "hgvsp_short"})
    .set_index("hgvsp_short", drop=True)
)
mut_sets = {}
for mut_set in mut_classes.columns:
    mut_sets[mut_set] = mut_classes.query(f"{mut_set} == True").index.to_list()

#########################################################################
# Setup
#########################################################################

# Where to find the GENIE data
genie_dir = os.path.join(os.environ.get("HOME"), "genie/Data/Original/genie-15.0/")

# Initialize the Genie object
g = gd.Genie(genie_dir)

#########################################################################
# Define what we are interested in
#########################################################################

## Gene of interest
gene_symbols = ["ERBB2"]

# We want to focus on MANE transcripts
universe = "mane"

# Cancer types we we are interested in
nsclc = "Non-Small Cell Lung Cancer"
crc = "Colorectal Cancer"
blca = "Bladder Cancer"
brca = "Breast Cancer"
cancer_types = [nsclc, crc, blca, brca]

# For NSCLC, we want to keep only LUAD, while for the other indications we
# don't differentiate by subtype. There are two possible approaches how to deal
# with that:
# (1) Work on subtype level, filter results for subtype "LUAD" for NSCLC, and
#     for subtype "all" for all other indications.
# (2) Identify all sample ids for indications of interest, that is, only LUAD
#     samples for NSCLC and all samples for other indications, and specify the
#     sample ids when calling functions of the genie package.
# Here we implement option (1) because the required code is shorter, although
# the intermediate memory requirements and runtime may be larger.

cancer_subtypes_to_keep = ["LUAD"]

#########################################################################
# Get unique mutations found in GENIE, filtered for colleague's mutation sets
#########################################################################

muts = (
    g.get_annotated_unique_mutations(gene_symbols, universe)
    .query("MANE_status=='MANE Select'")
    .loc[:, ["hgvsg", "hgvsp", "hgvsp_short"]]
)
hgvsgs = {}
for mut_set, hgvsp_shorts in mut_sets.items():
    df = muts.query("hgvsp_short.isin(@hgvsp_shorts)")
    hgvsgs[mut_set] = df.hgvsg.to_list()

#########################################################################
# Get gene level counts and frequencies
#########################################################################

idx = [(x, "all") for x in cancer_types if x != nsclc]
idx.append((nsclc, "LUAD"))
freqs_gl = {}
for mut_set in mut_sets.keys():
    freqs_gl[mut_set] = (
        g.get_gene_level_frequencies(
            gene_symbols=gene_symbols,
            hgvsgs=hgvsgs[mut_set],
            universe=universe,
            cancer_types=cancer_types,
            cancer_subtypes_to_keep=cancer_subtypes_to_keep,
            panel_coverage_threshold=0.8,
            impute=True,
            precision=3,
        )
        .droplevel(level="gene_symbol")
        .loc[idx]
    )

#########################################################################
# Write results to an Excel file
#########################################################################

with pd.ExcelWriter("HER2_GENIE_15.0_selected_cancers_gene_level_imputed.xlsx") as writer:
    for mut_set, freqs in freqs_gl.items():
        freqs.sort_index().to_excel(writer, sheet_name=mut_set, merge_cells=False)
