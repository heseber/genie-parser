#!/usr/bin/env python
import gzip
import hashlib
import os
import re

import pandas as pd

mut_file = "../Data/Original/genie-15.0/data_mutations_extended.txt"
vcf_file = "../Data/Original/genie-15.0/genie.vcf.gz"


class GenPosChecker:
    # This class is needed for adding missing leftmost bases of indels
    def __init__(self, gencode_dir=os.path.join(os.getenv("HOME"), "gencode", "chr")):
        self.chromosome_cache = dict()
        self.gencode_dir = gencode_dir

    def _find_chr_file(self, chromosome):
        c = re.sub("chr", "", str(chromosome))
        if c == "23":
            c = "X"
        if c == "MT":
            c = "M"
        pattern = f".*chr{c}.fa.gz"
        files = os.listdir(self.gencode_dir)
        file = [x for x in files if re.match(pattern, x)][0]
        file = os.path.join(self.gencode_dir, file)
        return file

    def _read_chr(self, file_name):
        with gzip.open(file_name, "rt") as file:
            # Read the line with the sequence name
            file.readline()
            # Read the sequence
            seq = "".join(x.strip() for x in file)
            return seq

    def get_chr_seq(self, chromosome):
        c = str(chromosome)
        if c not in self.chromosome_cache:
            chr_file = self._find_chr_file(c)
            self.chromosome_cache[c] = self._read_chr(chr_file)
        return self.chromosome_cache[c]

    def get_chr_sub_seq(self, chromosome, start, stop, first_pos=1):
        seq = self.get_chr_seq(chromosome)
        start = max(0, start - first_pos)
        stop = min(len(seq), stop - first_pos + 1)
        return seq[start:stop]

    def get_reverse_complement(self, seq):
        compl = {"A": "T", "T": "A", "C": "G", "G": "C"}
        seq = seq[::-1]
        seq = "".join([compl[b] for b in seq])
        return seq

    def get_strand(self, chromosome, pos, ref, prefix_len=0, first_pos=1):
        ref_rc = self.get_reverse_complement(ref)
        pos = pos - prefix_len - first_pos
        seq = self.get_chr_seq(chromosome)
        subseq = seq[pos : (pos + len(ref))]
        if subseq == ref:
            return "+"
        elif subseq == ref_rc:
            return "-"
        else:
            return None


def get_sha256(row):
    vid = ":".join(
        [
            row.Chromosome,
            str(row.Start_Position),
            row.Reference_Allele,
            row.Tumor_Seq_Allele2,
        ]
    )
    return hashlib.sha256(vid.encode("UTF-8")).hexdigest()


def long_chromosome_names(series):
    s = "chr" + series.map(str)
    s = s.str.replace("chrMT", "chrM")
    return s


gpc = GenPosChecker()

df = pd.read_table(mut_file, low_memory=False)
df = df[["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]]
df = df.drop_duplicates()
df["ID"] = df.apply(get_sha256, axis=1)


# Fix indels with missing leftmost common nucleotide
df["Reference_Allele"] = df["Reference_Allele"].replace("-", "")
df["Tumor_Seq_Allele2"] = df["Tumor_Seq_Allele2"].replace("-", "")
for idx, row in df.iterrows():
    if row.Reference_Allele and row.Tumor_Seq_Allele2:
        continue
    # As it seems, if no reference allele was specified for an insertion
    # (i.e., Reference_Allele == "-"), then the Start_Position refers to
    # the nucleotide left of the insertion, so this is the position _after_
    # which the insertion is added. The Start_Position and End_Position then
    # are the positions between which the sequence is inserted. To get a
    # standard VCF, the reference allele needs to be replaced with the base at
    # Start_Position, no shift of the Start_Position is needed in that case.
    # For deletions, the position left of the deletion needs to be added, so
    # the Start_Position needs to be decremented and the base at that position
    # needs to be added to the reference and alternative alleles.
    offset = 0
    if not row.Tumor_Seq_Allele2:  # deletion
        offset = 1
    leading_base = gpc.get_chr_sub_seq(
        row.Chromosome, row.Start_Position - offset, row.Start_Position - offset
    )
    df.loc[idx, "Reference_Allele"] = leading_base + df.loc[idx, "Reference_Allele"]
    df.loc[idx, "Tumor_Seq_Allele2"] = leading_base + df.loc[idx, "Tumor_Seq_Allele2"]
    df.loc[idx, "Start_Position"] -= offset


vcf = pd.DataFrame(
    {
        "#CHROM": long_chromosome_names(df.Chromosome),
        "POS": df.Start_Position,
        "ID": df.ID,
        "REF": df.Reference_Allele,
        "ALT": df.Tumor_Seq_Allele2,
        "QUAL": ".",
        "FILTER": "PASS",
        "INFO": ".",
        "FORMAT": "GT:VF",
        "GENIE": "1/1:1.0",
    }
).sort_values(["#CHROM", "POS"])

# Remove lines where REF and ALT are identical (in genie 15.0, there are 4
# such entries)
vcf = vcf[vcf.REF != vcf.ALT]

os.makedirs(os.path.dirname(vcf_file), exist_ok=True)
with gzip.open(vcf_file, "wt") as f:
    f.writelines(
        [
            "##fileformat=VCFv4.3\n",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            '##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant Frequency">\n',
        ]
    )
    vcf.to_csv(f, sep="\t", index=False)
