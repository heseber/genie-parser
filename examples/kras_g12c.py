#!/usr/bin/env python

# Runtime of this script are about 1min 12s

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
gene_symbols = ["KRAS"]

# We want to focus on MANE transcripts
universe = "mane"

# Cancer types we we are interested in
nsclc = "Non-Small Cell Lung Cancer"
crc = "Colorectal Cancer"
cancer_types = [nsclc, crc]

# Cancer subtypes not to summerize in the "other" category
cancer_subtypes_to_keep = ["LUAD", "LUSC", "COAD", "COADREAD", "READ"]

# We want to have frequencies for the different races, and we want to keep
# "White", "Black" and "Asian", while summarizing all other races in "other".
extra_group_columns = {"PRIMARY_RACE": ["White", "Black", "Asian"]}

# Get the mutation leading to G12C
muts = (
    g.get_annotated_unique_mutations(gene_symbols, universe)
    .query("MANE_status=='MANE Select'")
    .query("hgvsp.notna()")
    .query("hgvsp_short=='G12C'")
    .loc[:, ["hgvsg", "hgvsp", "hgvsp_short"]]
)
hgvsgs = muts.hgvsg.to_list()
hgvsps = muts.hgvsp.to_list()

################################################################################
# Get the mutation frequencies at protein sequence level
################################################################################

freqs = g.get_amino_acid_level_frequencies(
    gene_symbols=gene_symbols,
    hgvsgs=hgvsgs,
    hgvsps=hgvsps,
    cancer_types=cancer_types,
    cancer_subtypes_to_keep=cancer_subtypes_to_keep,
    extra_group_columns=extra_group_columns,
    universe=universe,
    precision=1,
).droplevel("hgvsp")


################################################################################
# Write Excel file with results
################################################################################


# For the specified cancer type and race, get a sheet with patient counts for
# WT, MUT and N, and the allele frequencies, including confidence intervals.
def get_sheet(freqs, cancer_type, race):
    df = (
        freqs.query("CANCER_TYPE==@cancer_type")
        .query("PRIMARY_RACE==@race")
        .droplevel("CANCER_TYPE")
        .droplevel("PRIMARY_RACE")
        .reset_index()
        .rename(columns={"ONCOTREE_CODE": "subtype"})
    )
    return df


# For the specified cancer type, get a summary sheet with formatted allele
# frequencies, including confidence intervals, with cancer subtypes in rows and
# race in columns.
def get_summary_sheet(freqs, cancer_type):
    df = (
        freqs.query("CANCER_TYPE==@cancer_type")["AF_FORMATTED"]
        .unstack("PRIMARY_RACE")[["White", "Black", "Asian", "other", "all"]]
        .droplevel("CANCER_TYPE")
        .reset_index()
        .rename(columns={"ONCOTREE_CODE": "subtype"})
    )
    return df


# Finally, write the excel file
races = ["White", "Black", "Asian", "other", "all"]
with pd.ExcelWriter("KRAS_G12C_GENIE_15.0.xlsx", engine="xlsxwriter") as writer:
    get_summary_sheet(freqs, nsclc).to_excel(
        writer, sheet_name="NSCLC Summary", index=False
    )
    get_summary_sheet(freqs, crc).to_excel(
        writer, sheet_name="CRC Summary", index=False
    )
    for race in races:
        get_sheet(freqs, nsclc, race).to_excel(
            writer, sheet_name="NSCLC " + race, index=False
        )
    for race in races:
        get_sheet(freqs, crc, race).to_excel(
            writer, sheet_name="CRC " + race, index=False
        )
