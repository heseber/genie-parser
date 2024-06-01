#!/bin/bash

DIR="../Data/Original/consortium/Main_GENIE_cBioPortal_Releases/16.2-consortium"

ENV="genie"

__conda_setup="$('$HOME/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
        . "$HOME/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="$HOME/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "$HOME/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "$HOME/miniforge3/etc/profile.d/mamba.sh"
else
    echo "Please install mamba first."
    echo "See https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html"
    exit 1
fi

if ! mamba env list | grep -qe "^$ENV "; then
    echo "The environment '$ENV' does not exists"
    echo "Please create this environment!"
    exit 1
fi

mamba activate $ENV

if ! type bcftools >/dev/null 2>&1; then
    echo "Please install bcftools in first."
    echo "mamba activate nirvana"
    echo "mamba install bcftools"
    exit 1
fi

FA=$HOME/gencode/GRCh37.primary_assembly.genome.fa.bgz
if ! test -f $FA; then
    echo "'$FA' is missing, please dowload and install!"
    exit 1
fi

FAI=$HOME/gencode/GRCh37.primary_assembly.genome.fa.bgz.fai
if ! test -f $FAI; then
    echo "'$FAI' is missing, please create the FASTA index!"
    exit 1
fi

cd $DIR
zcat genie.vcf.gz | bgzip > tmp.vcf.bgz
bcftools reheader --fai $FAI tmp.vcf.bgz -o tmp2.vcf.bgz
bcftools norm -f $FA -c w -o genie_normalized.vcf.bgz -Oz tmp2.vcf.bgz
rm tmp.vcf.bgz tmp2.vcf.bgz
