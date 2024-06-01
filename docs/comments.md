# Comments on the GENIE data

## Mutations are not normalized

The mutations contained in GENIE are not always normalized (left-shifted).
Insertions and deletions do not always include the leading common nucleotide.
Therefore, one and the same genomic alteration may be reported differently by
different assays. When mutations are aggregated across different panels, the
specification of genomic alterations must be normalized.

This is done by creating a VCF file from the GENIE mutation data, normalizing
and re-annotating the VCF file (here we use *Illumina Connected Annotations*).

## GENIE uses ENSEMBL transcripts

Although the annotation of mutations provided by GENIE contains also RefSeq
transcripts, the RefSeq annotation is incomplete. Many mutations have only an
ENSEMBL transcripts but no RefSeq transcript assigned. On the other hand, there
is no alteration that has only an assigned RefSeq transcript but no ENSEMBL
transcript.

Based on this observation, there are several different options for creating a
unified consistent annotation across panels:

### Option 1: Keep ENSEMBL transcripts reported by GENIE

Consistently re-annotate all genomic variants, extract annotations for the
ENSEMBL transcripts contained in GENIE, and map these updated annotations to the
genomic changes reported by GENIE, using the ENSEMBL transcript ID as a mapping
link.

This option keeps the annotation as close a possible to the original GENIE
annotations and just ensures that we get consistent annotations across all
assays.

However, this option does not have any focus on MANE Select transcripts.

### Option 2: Focus on MANE transcripts

Clinical reporting of mutations should preferrably be based on MANE transcripts.
GENIE is based on GRCh37, while MANE is defined for GRCh38. The NCBI has
remapped the MANE transcript RefSeq sequences back to the GRCh37 genome. Ensembl
did not do this. Although [there is a mapping table](
http://tark.ensembl.org/web/mane_GRCh37_list/ ) from MANE transcripts to Ensembl
transcripts of GRCh37, this does not cover all MANE transcripts and the GRCh37
transcripts are also not always identical to the MANE transcripts. Therefore,
when focusing on MANE transcripts while working with the GRCh37 genome (which
GENIE is based on), using RefSeq transcripts rather than Ensembl transcripts is
preferrable.

Focusing on MANE requires some additional decisions. If a genomic variant
affects a MANE transcript, the consequence on that MANE transcript can be
reported. Even then, it needs to be decided how to handle situations where both
*MANE Select* and *MANE Plus Clincial* transcripts are affected.

If a genomic variant does not affect a MANE transcript, there is no clear and
obvious procedure as to which transcript to report. One could decide to report
only effects on MANE transcripts and omit all genomic variants that do not
affect a MANE transcript. Or one could decide to omit non-MANE transcripts for
genes for which a MANE transcript is defined, even if the MANE transcript itself
is not affected by a genomic variant, but to report none-MANE transcripts for
genes for which no MANE transcript is defined. For such genes it still must be
decided which of the non-MANE transcripts to report.

So there are several different options when focusing on MANE. Unless the report
is MANE-only, i.e. filters out any non-MANE transcript, it depends very much on
the use case which transcripts to keep.

### Conclusion

The `genie` package provides two types of annotations, one is based on the
ENSEMBL transcripts originally used by Genie, the other one is MANE-only and
uses RefSeq. The `universe` argument of many `genie` functions can be used to
select the annotation type.

If an enhanced "MANE + Non-MANE" version is required, the file
`genie.annot.tsv.gz` can be used as input for a customized filtering process
outside of the `genie` package. This file was generated when creating the
auxiliary files for a new GENIE release and contains annotations for all
transcripts affected by a genomic variant.  



