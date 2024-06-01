# Setting up the environment

## Introduction

Some preparations are required before downloading the GENIE data and  creating
additional auxiliary files needed by the `genie` package. Please follow the
instructions below to set up your environment.

## Install mamba

If you haven't done so already anyway, please install `mamba`, the next
generation replacement for `conda`. One of the advantages of `mamba` is that
it is much faster than `conda`.

Please see the [Mamba Installation Instructions](
    https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
) on how to install `mamba`.

## Create a virtual environment for genie

!!! example "Create and activate a virtual environment"

    ```sh
    mamba create -n genie
    mamba activate genie
    ```

## Install genie package

!!! example "Install genie package"

    ```sh
    mamba install python pandas sinfo
    pip install "git+https://github.com/Bayer-Group/genie.parser.git#subdirectory=genie"
    ```

## Install packages and programs


!!! example "Install packages"

    ```sh title="Installation"
    mamba install -c bioconda bcftools samtools
    pip install synapseclient
    ```

## Install Illumina Annotator and nirvanaparser module

For the normalization and annotation of mutations from GENIE, a working setup of
the Illumina Annotator (a.k.a. Nirvana) and the in-house developed
 `nirvanaparser` model is needed. Setting up this environment is beyond the
 scope of this document, please see the [documentation of `nirvanaparser`](
 https://github.com/Bayer-Group/BDS-MutationAnnotationFiltering/tree/main/nirvana_parser
 ) for instructions.