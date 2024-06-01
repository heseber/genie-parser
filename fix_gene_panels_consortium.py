#!/usr/bin/env python
# coding: utf-8

import os

import pandas as pd

genie_dir = "../Data/Original/consortium/Main_GENIE_cBioPortal_Releases/16.2-consortium"


class Samples:
    def __init__(self, genie_dir):
        self.sample_info = pd.read_table(
            os.path.join(genie_dir, "data_clinical_sample.txt"),
            index_col="SAMPLE_ID",
            comment="#",
            low_memory=False,
        )

    def get_sample_info(self, sample_id: str) -> dict:
        return self.sample_info.loc[sample_id].to_dict()


samples = Samples(genie_dir)
mut = pd.read_table(
    os.path.join(genie_dir, "data_mutations_extended.txt"), low_memory=False
)

wxs_samples = set(samples.sample_info.query('SEQ_ASSAY_ID=="PROV-TRISEQ-V2"').index)
wxs_genes = set(mut.loc[mut.Tumor_Sample_Barcode.isin(wxs_samples), "Hugo_Symbol"])
fname = os.path.join(genie_dir, "gene_panels", "data_gene_panel_PROV-TRISEQ-V2.txt")
with open(fname, "wt") as fh:
    fh.write("stable_id: PROV-TRISEQ-V2\n")
    fh.write(f"description: PROV-TRISEQ-V2, Number of Genes - {len(wxs_genes)}\n")
    fh.write("\t".join(["gene_list:"] + list(wxs_genes)) + "\n")
