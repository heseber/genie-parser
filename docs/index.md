# Welcome to the GENIE Python package



## Introduction

The GENIE project provides mutation data for a large number of patients and
cancer indications. However, many different panels are used for generating the
data, and these panels have different coverage of tested genomic ranges.
Therefore, to get the mutation profile of genes across samples it needs to be
checked whether a particular genomic region is tested by a panel. If a sample is
tested with a panel not covering a region of interest, the mutational status of
this sample is unknown.

The `genie` package provides means to get mutation profiles of genes for all
samples where the gene has actually been tested. Such profiles can then also be
used to estimate the frequency of mutations in a disease indication.

Please see also [some comments on the GENIE data](comments.md) for observations
that I made while working on this package.

