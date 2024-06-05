#!/usr/bin/env python
# Author: Henrik Seidel <henrik.seidel@bayer.com>
# Date: 2024-03-20


import os

import genie.genie as gd
import pandas as pd

################################################################################
# Setup
################################################################################

# Where to find the GENIE data
genie_dir = os.path.join(os.environ.get("HOME"), "genie/Data/Original/genie-15.0/")

# Initialize the Genie object
g = gd.Genie(genie_dir)

################################################################################
# Define what we are interested in
################################################################################

# Gene of interest
gene_symbols = ["ERBB2"]

# We want to focus on MANE transcripts
universe = "mane"

# We focus on the top 5 cancer types with the most frequent mutations
cancer_types = [
    "Bladder Cancer",
    "Small Bowel Cancer",
    "Cervical Cancer",
    "Non-Small Cell Lung Cancer",
    "Esophagogastric Cancer",
]

# These are the cancer subtypes we want not to be summarized in the "other" category
cancer_subtypes_to_keep = [
    "BLCA",
    "UTUC",
    "BLAD",
    "UCU",
    "DA",
    "SBC",
    "SIC",
    "CESC",
    "ECAD",
    "CEAD",
    "LUAD",
    "LUSC",
    "ESCA",
    "STAD",
    "GEJ",
    "EGC",
    "ESCC",
    "SPDAC",
    "SSRCC",
    "DSTAD",
    "EPDCA",
]

# We want to get the frequencies by race, and any race except "White",
# "Black" and "Asian" will be summarized in "other".
extra_group_columns = {"PRIMARY_RACE": ["White", "Black", "Asian"]}

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

################################################################################
# Get the frequencies of amino acid variants
################################################################################

freqs_aa = g.get_amino_acid_level_frequencies(
    gene_symbols=gene_symbols,
    hgvsgs=hgvsgs,
    hgvsps=hgvsps,
    cancer_types=cancer_types,
    cancer_subtypes_to_keep=cancer_subtypes_to_keep,
    extra_group_columns=extra_group_columns,
    universe=universe,
    precision=2,
)
# Add hgvsp_short to index, which is okay because we have focused on a single
# amino acid sequence for HER2
freqs_aa = (
    freqs_aa.reset_index()
    .merge(muts[["hgvsp", "hgvsp_short"]].drop_duplicates(), on="hgvsp")
    .set_index(["hgvsp_short"] + freqs_aa.index.names, drop=True)
)

################################################################################
# Print and save results
################################################################################

# Initialize list of Excel file sheets
sep = "\n" + "-" * 80 + "\n\n"
sheets = {}

# All top 5 mutations

print(sep + "MUTATIONS FOUND IN AT LEAST 10 PATIENTS PER INDICATION:")
df = (
    freqs_aa.query("MUT>10")
    .query('CANCER_TYPE!="all"')
    .query('ONCOTREE_CODE=="all"')
    .query('PRIMARY_RACE=="all"')
    .sort_values("AF_PERC", ascending=False)
)
print(df)
sheets["AF for MUT>=10"] = df.copy()

# NSCLC

print(sep + "\nTOP 10 MUTATIONS IN NSCLC:")
df = (
    freqs_aa.query('CANCER_TYPE=="Non-Small Cell Lung Cancer"')
    .query('PRIMARY_RACE=="all"')
    .query('ONCOTREE_CODE=="all"')
    .sort_values("AF_PERC", ascending=False)
    .head(10)
)
print(df)
nsclc_top5_muts = df.head(5).index.get_level_values("hgvsp_short").to_list()
sheets["NSCLC top 10"] = df.copy()


print(sep + "TOP 5 MUTATIONS IN NSCLC BY RACE AND CANCER SUBTYPE:\n")
df = (
    freqs_aa.query('CANCER_TYPE=="Non-Small Cell Lung Cancer"')
    .query("hgvsp_short.isin(@nsclc_top5_muts)")
    .loc[:, "AF_FORMATTED"]
    .unstack("PRIMARY_RACE")[["White", "Black", "Asian", "other", "all"]]
)
print(df)
sheets["NSCLC top 5 by subtype+race"] = df.copy()

# Bladder cancer

print(sep + "TOP 10 MUTATIONS IN BLADDER CANCER:")
df = (
    freqs_aa.query('CANCER_TYPE=="Bladder Cancer"')
    .query('PRIMARY_RACE=="all"')
    .query('ONCOTREE_CODE=="all"')
    .sort_values("AF_PERC", ascending=False)
    .head(10)
)
print(df)
blca_top5_muts = df.head(5).index.get_level_values("hgvsp_short").to_list()
sheets["Bladder top 10"] = df.copy()

print(sep + "TOP 5 MUTATIONS IN BLADDER CANCER BY RACE AND CANCER SUBTYPE:\n")
df = (
    freqs_aa.query('CANCER_TYPE=="Bladder Cancer"')
    .query("hgvsp_short.isin(@blca_top5_muts)")
    .loc[:, "AF_FORMATTED"]
    .unstack("PRIMARY_RACE")[["White", "Black", "Asian", "other", "all"]]
)
print(df)
sheets["Bladder top 5 by subtype+race"] = df.copy()

# Write excel file

with pd.ExcelWriter("HER2_GENIE_15.0_top_cancers.xlsx", engine="xlsxwriter") as writer:
    for sheet_name, df in sheets.items():
        df.to_excel(writer, sheet_name=sheet_name, merge_cells=False)
