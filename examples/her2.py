#!/usr/bin/env python

##### THIS SCRIPT NEEDS >32GB <64GB of RAM AND RUNS FOR ABOUT 30 MINUTES #####

import os

import genie.genie as gd
import pandas as pd

###################################
# Setup
###################################

# Where to find the GENIE data
genie_dir = os.path.join(os.environ.get("HOME"), "genie/Data/Original/genie-15.0/")

# Initialize the Genie object
g = gd.Genie(genie_dir)

###################################
# Define what we are interested in
###################################

## Gene of interest
gene_symbols = ["ERBB2"]

# We want to focus on MANE transcripts
universe = "mane"

# We want to keep only protein modifying mutations for MANE Select transcripts
muts = (
    g.get_annotated_unique_mutations(gene_symbols, universe)
    .query("MANE_status=='MANE Select'")
    .query("consequence_priority <= 6")  # only protein sequence modifying mutations
    .query("consequence_priority > 4")  # no potentially deleterious mutations
    .loc[:, ["hgvsg", "hgvsp", "hgvsp_short"]]
)
hgvsgs = muts.hgvsg.to_list()
hgvsps = muts.hgvsp.to_list()

###############################################
# Get amino acid level counts and frequencies
###############################################

freqs_aa = (
    g.get_amino_acid_level_frequencies(
        gene_symbols=gene_symbols,
        hgvsgs=hgvsgs,
        hgvsps=hgvsps,
        universe=universe,
        panel_coverage_threshold=0.8,
        impute=False,
        precision=3,
    )
    .reset_index()
    .merge(
        g.mutations.annot[universe][["hgvsp", "hgvsp_short"]].drop_duplicates(),
        on="hgvsp",
    )
    .set_index(["hgvsp_short", "CANCER_TYPE"], drop=True)
    .drop(columns="hgvsp")
)

#########################################
# Get gene level counts and frequencies
#########################################

freqs_gl = g.get_gene_level_frequencies(
    gene_symbols=gene_symbols,
    hgvsgs=hgvsgs,
    universe=universe,
    panel_coverage_threshold=0.8,
    impute=False,
    precision=3,
)

###################################
# Write Excel file
###################################

with pd.ExcelWriter("HER2_GENIE_15.0_all_cancers.xlsx") as writer:
    freqs_aa.to_excel(writer, sheet_name="FREQS_AA", merge_cells=False)
    freqs_aa.query('CANCER_TYPE=="all"').sort_values(
        "AF_PERC_CI_LOWER", ascending=False
    ).to_excel(writer, sheet_name="FREQS_AA global", merge_cells=False)
    freqs_aa.query('CANCER_TYPE!="all"').query("N>=100").sort_values(
        "AF_PERC_CI_LOWER", ascending=False
    ).to_excel(writer, sheet_name="FREQS_AA by cancer for N>=100", merge_cells=False)
    freqs_gl.sort_values("AF_PERC_CI_LOWER", ascending=False).to_excel(
        writer, sheet_name="FREQS_Gene_Level", merge_cells=False
    )
    freqs_gl.sort_values("AF_PERC_CI_LOWER", ascending=False).query("N>=100").to_excel(
        writer, sheet_name="FREQS_Gene_Level for N>=100", merge_cells=False
    )
