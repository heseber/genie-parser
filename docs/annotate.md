# Annotating mutations

The panels included in GENIE do not consistently use left-shifted specifications
of genomic variants. As a consequence, two samples that were tested with
different panels may carry the same mutation, but the specification of the
mutation (chromosome, position, reference allele, alternative allele) can be
different. Another issue is that some panels report a "-" in the reference
allele column for insertions or for the alternative allele for deletions. This
is not consistent with the VCF file format standard and breaks downstream tools. 

For these reasons, the specification of mutations must be normalized by adding
the missing leftmost nucleotide of indels, and left-shifting mutations,
resulting in HGVSG designations of the normalized variants.

This is a multi-step process and needs to be done only once for every new GENIE
release.

## Step 1: Prepare the reference genome

For adding the missing nucleotide left of inserations or deletions, the genomic
sequence is required.

### Download reference sequence

!!! example "download_reference.sh"

    ```sh
    GENCODE_DIR=$HOME/gencode
    mkdir -p $GENCODE_DIR
    cd $GENCODE_DIR
    FILE="GRCh37.primary_assembly.genome.fa.gz"
    URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/\
    latest_release/GRCh37_mapping/$FILE"
    curl $URL > $HOME/$FILE
    ```

### Bgzip and index reference sequence

!!! example "reformat_genome.sh"

    ```sh
    zcat GRCh37.primary_assembly.genome.fa.gz \
    | bgzip > GRCh37.primary_assembly.genome.fa.bgz
    samtools faidx GRCh37.primary_assembly.genome.fa.bgz
    ```

### Create chromosome specific FASTA files

A subdirectory with one FASTA file for each chromosome needs to be provided. The
genome FASTA file can be split into chromosome specific files using the Python
script `extract_sequences.py` which is part of the GitHub repository of the
Bayer GENIE project and also included here:

??? note "extract_sequences.py"
    ```python
    # Author: Henrik Seidel
    # Date: 2018-07-03
    # This script extracts all sequences from a fasta file into single sequences

    import sys
    import os.path
    import gzip
    import re

    if len(sys.argv) != 3:
        sys.exit("ERROR: you have to specify the name of the fasta file " +
        "and the name of the output directory")

    file_name = sys.argv[1]
    out_dir = sys.argv[2]

    # Assert we have a filename as first argument
    if not(os.path.isfile(file_name)):
        sys.exit("ERROR: No such file: \"" + file_name + "\"")

    # Assert the name of the output directory is okay
    if not(os.path.exists(out_dir)):
        os.mkdir(out_dir)
    elif not(os.path.isdir(out_dir)):
        sys.exit("ERROR: \"" + out_dir + "\" exists but is not a directory")

    # Choose function for opening file, depending on compression status of file
    if (file_name.endswith(".gz") or file_name.endswith(".bgz")):
        this_open = gzip.open
        file = gzip.open(file_name, mode = "rt")
    else:
        this_open = open

    n = 0
    with this_open(file_name, mode = "rt") as file:
        out_file = None
        try:
            for line in file:
                if line.startswith(">"):
                    if out_file is not None:
                        out_file.close()
                    seq_name = re.match("^>([^ ]+)", line)[1]
                    n += 1
                    out_file_name = '{:s}/{:03d}_{:s}.fa'.format(out_dir, n, seq_name)
                    print('Writing {:s} to {:s}'.format(seq_name, out_file_name))
                    out_file = open(out_file_name, mode = "wt")
                out_file.write(line)
        finally:
            out_file.close()
    ```

!!! example "split_genome.sh"

    ```sh
    mamba activate genie
    python extract_sequences.py GRCh37.primary_assembly.genome.fa.bgz chr
    bzgip chr/*.fa
    ```

Required directory structure:

??? info "Directory structure of genome sequence files"
    ```
    $HOME/gencode
    ├── GRCh37.primary_assembly.genome.fa.bgz
    ├── GRCh37.primary_assembly.genome.fa.bgz.fai
    ├── GRCh37.primary_assembly.genome.fa.bgz.gzi
    ├── GRCh37.primary_assembly.genome.fa.gz
    └── chr
        ├── 001_chr1.fa.gz
        ├── 002_chr10.fa.gz
        ├── 003_chr11.fa.gz
        ├── 004_chr12.fa.gz
        ├── 005_chr13.fa.gz
        ├── 006_chr14.fa.gz
        ├── 007_chr15.fa.gz
        ├── 008_chr16.fa.gz
        ├── 009_chr17.fa.gz
        ├── 010_chr18.fa.gz
        ├── 011_chr19.fa.gz
        ├── 012_chr2.fa.gz
        ├── 013_chr20.fa.gz
        ├── 014_chr21.fa.gz
        ├── 015_chr22.fa.gz
        ├── 016_chr3.fa.gz
        ├── 017_chr4.fa.gz
        ├── 018_chr5.fa.gz
        ├── 019_chr6.fa.gz
        ├── 020_chr7.fa.gz
        ├── 021_chr8.fa.gz
        ├── 022_chr9.fa.gz
        ├── 023_chrM.fa.gz
        ├── 024_chrX.fa.gz
        └── 025_chrY.fa.gz
    ```

## Step 2: Generate a VCF file of all GENIE mutations

A VCF file is required for annotating mutations. All mutations from the GENIE
file `data_mutations_extended.txt` will be extracted and transformed to a VCF
file by the script [`create_vcf.py`](
https://github.com/Bayer-Group/genie.parser/blob/main/create_vcf.py ) which is
part of the GitHub genie project.

!!! warning

    Before calling `create_vcf.py` you need to edit the first lines where
    directory locations are defined.

!!! example "create_vcf.sh"

    ```sh
    # Please be patient, this takes about 6 minutes
    mamba activate genie
    python create_vcf.py
    ```

## Step 3: Normalize the VCF file

Normalizing the VCF file before annotating it with *Illumina Connected
Annotations* (previously known as *Nirvana*) is currently necessary because
otherwise the mapping between samples and mutations is lost. This will change
with forthcoming releases of *Illumina Connected Annotations* which will
transfer the variant identifiers from the VCF ID column to the JSON output file
of the annotator.

!!! warning

    Before calling `normalize_vcf.sh` you need to edit the first lines where
    the directory location and the environment name are defined.

!!! example "normalize_vcf.sh"

    ```sh
    ./normalize_vcf.sh
    ```

## Step 4: Annotate mutations from the normalized VCF file

The VCF file needs to be annotated with *Illumina Connected Annotations* (previously known as *Nirvana*). See the [documentation of the `icaparser` package](
https://github.com/Bayer-Group/ica-parser) on how to set up and use the annotator. The `genie` GitHub project contains two
scripts for running the annotator - `run_ica.sh` and
`run_ica_docker.sh`. Use one of these two scripts for annotating the VCF
file, depending on how you installed the annotator.

!!! example "Example for running the annotator"

    ```sh
    # This takes about 5 minutes
    ./run_ica_docker.sh ../Data/Original/genie-15.0/genie_normalized.vcf.bgz
    mv ../Data/Original/genie-15.0/genie_normalized.json* ../Results
    ```

## Step 5: Create auxiliary files from annotations

The results of the annotator will be used to create two auxiliary files:

1. `genie.vcf_id_to_hgvsg.tsv.gz`: This file contains a mapping between variant
   identifiers and the normalized hgvsg specification of mutations. 
2. `genie.annot.tsv.gz`: This file contains annotation of all variants from the
   VCF file.

These files are created by the **Jupyter notebook
[`genie_create_aux_files.ipynb`](
https://github.com/Bayer-Group/genie-parser/genie_create_aux_files.ipynb
)** which is part of of the genie GitHub project. Please open this file in a
Jupyter server, adjust directory names as needed, and run it. 

!!! warning

    Running the notebook will take about one hour, and requires a lot of
    memory. Make sure you use an instance with 64 GB RAM or more, and no memory
    intensive processes running in parallel.

!!! note

    All required files are saved in the `Data/Original` directory subtree,
    although some of the files are derived files. This is not according to the
    standard guidelines. The rationale for this procedure is that in the end
    we want to have a single directory with all data files and auxiliary files
    needed by the `genie` Python module.

## Step 6: Precompute which mutations are tested by which panel

GENIE uses different panels (assays), and for determining mutation profiles and
mutation frequencies it is important to know whether or not a panel tests a
particular mutation. Determining this on the fly is extremely time consuming.
Therefore, this needs to be precomputed for each new GENIE release.

To be more specific, a big matrix will be precomputed with panels versus
mutations, where mutations are all mutations detected in at least one sample.
This matrix is then saved as a Parquet file `mutation_tested_by_panel.parquet`.

This Python script makes use of the `genie` package. Because the `Programs`
directory where the script `create_panel_tested_positions_cache.py` is located
also contains the sources of the `genie` package in the subdirectory `genie`,
it is necessary to change the working directory to somewhere else, like the
home directory, when calling the script, to make sure that it picks up the
package as it was installed by the previous Python environment setup step. If
you fail to do so, you will get an error message concering the `genie` package.

!!! example "Precompute mutation_tested_by_panel.parquet"

    ```sh
    # This takes about a day
    cd $HOME
    mamba activate genie
    SCRIPT=$HOME/genie/Programs/create_panel_tested_positions_cache.py
    python $SCRIPT
    ```

!!! warning

    This calculation takes about a full day. Make sure the S@S EZ instance has
    the stop scheduler turned off during that time.