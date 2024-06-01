#!/bin/bash
NDIR=$HOME/ICA
GENOME=GRCh37
VCF=$1
BASE=$(basename $VCF .vcf)
echo "*********************************************************"
echo "*  Annotating $VCF"
echo "*********************************************************"
dotnet $NDIR/Nirvana.dll \
-c   $NDIR/Data/Cache/ \
--sd $NDIR/Data/SupplementaryAnnotation/$GENOME \
-r   $NDIR/Data/References/Homo_sapiens.$GENOME.Nirvana.dat \
-i   $VCF \
-o   $BASE
