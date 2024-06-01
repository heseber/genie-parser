#!/bin/bash
N_DATA_DIR=$HOME/ICA
GENOME=GRCh37
VCF=$(realpath $1)
BASE=${VCF%%.vcf.gz}
BASE=${BASE%%.vcf.bgz}
echo "*********************************************************"
echo "*  Annotating $(basename $VCF)"
echo "*********************************************************"
docker run -it --rm -v $HOME:$HOME illumina-connected-annotations:3.22.0 Nirvana \
-c   $N_DATA_DIR/Cache/ \
--sd $N_DATA_DIR/SupplementaryAnnotation/$GENOME \
-r   $N_DATA_DIR/References/Homo_sapiens.$GENOME.Nirvana.dat \
-i   $VCF \
-o   $BASE
