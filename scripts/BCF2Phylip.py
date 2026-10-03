#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# BCF2Phylip.py
#
# Derived from vcf2phylip v2.9 by Edgardo M. Ortiz, with credits to
# Juan D. Palacio-Mejía.
#
# Modified 2026 by Adhityo Wicaksono / BCF2Phylip prototype project:
#   - direct BCF input through pysam/HTSlib
#   - FORMAT-aware GT access rather than assuming GT is first
#   - integer allele-index handling
#   - VCF2Phylip-compatible nucleotide/IUPAC semantics
#
# This modified work is distributed under the GNU General Public License
# version 3 (GPL-3.0), consistent with the upstream project.
#
# Upstream: https://github.com/edgardomortiz/vcf2phylip

"""
Convert a collection of SNPs in BCF format directly into PHYLIP, FASTA,
NEXUS, or binary NEXUS matrices for phylogenetic analysis.

This prototype intentionally follows the biological/output semantics of
vcf2phylip v2.9 while replacing text-VCF parsing with direct BCF access
through pysam/HTSlib.

Any ploidy is allowed for nucleotide matrices. Binary NEXUS is produced
only for diploid, biallelic SNP genotypes compatible with SNAPP encoding.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

__author__ = "Edgardo M. Ortiz; modified by Adhityo Wicaksono"
__credits__ = "Juan D. Palacio-Mejía"
__version__ = "0.2.0-prototype"
__date__ = "2026-10-04"
__upstream_version__ = "vcf2phylip 2.9"

# Kept exactly from vcf2phylip v2.9 for output compatibility.
# '*' is a deletion. Lowercase consensus is used when an N or * is part
# of an otherwise interpretable genotype.
AMBIG = {
    "A":"A", "C":"C", "G":"G", "N":"N", "T":"T",
    "*A":"a", "*C":"c", "*G":"g", "*N":"n", "*T":"t",
    "AC":"M", "AG":"R", "AN":"a", "AT":"W", "CG":"S",
    "CN":"c", "CT":"Y", "GN":"g", "GT":"K", "NT":"t",
    "*AC":"m", "*AG":"r", "*AN":"a", "*AT":"w", "*CG":"s",
    "*CN":"c", "*CT":"y", "*GN":"g", "*GT":"k", "*NT":"t",
    "ACG":"V", "ACN":"m", "ACT":"H", "AGN":"r", "AGT":"D",
    "ANT":"w", "CGN":"s", "CGT":"B", "CNT":"y", "GNT":"k",
    "*ACG":"v", "*ACN":"m", "*ACT":"h", "*AGN":"r", "*AGT":"d",
    "*ANT":"w", "*CGN":"s", "*CGT":"b", "*CNT":"y", "*GNT":"k",
    "ACGN":"v", "ACGT":"N", "ACNT":"h", "AGNT":"d", "CGNT":"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*":"-", "*ACGNT":"N",
}

GEN_BIN_TUPLES = {
    (None, None): "?",
    (0, 0): "0",
    (0, 1): "1",
    (1, 0): "1",
    (1, 1): "2",
}


def require_pysam():
    if pysam is None:
        print(
            "\nBCF2Phylip requires pysam/HTSlib.\n"
            "Install with:\n"
            "  python -m pip install pysam\n"
            "or:\n"
            "  conda install -c bioconda pysam\n"
        )
        sys.exit(1)


def extract_sample_names(variant_file):
    """Return BCF sample names in header order."""
    return [str(name).replace("./", "") for name in variant_file.header.samples]


def _normalize_allele(allele, ref):
    if allele is None:
        return "N"
    allele = str(allele).replace("-", "*").upper()
    if allele == "<NON_REF>":
        return ref
    return allele


def is_snp(record):
    """
    Match vcf2phylip's SNP rule:
    REF must have length 1, and each ALT (after <NON_REF> -> REF) must
    represent one character. Multiallelic SNPs are allowed.
    """
    ref = _normalize_allele(record.ref, "N")
    alts = record.alts or ()
    if len(ref) != 1 or not alts:
        return False
    normalized_alts = [_normalize_allele(a, ref) for a in alts]
    return all(len(a) == 1 for a in normalized_alts)


def num_genotypes(record, sample_names):
    """
    Match vcf2phylip v2.9 missingness semantics.

    In the text implementation a sample is considered missing when its sample
    field starts with '.', effectively making ./1 missing but 0/. present when
    GT is first. Here we reproduce that rule using the first GT allele.
    """
    missing = 0
    for name in sample_names:
        gt = record.samples[name].get("GT")
        if not gt or gt[0] is None:
            missing += 1
    return len(sample_names) - missing


def _allele_table(record):
    ref = _normalize_allele(record.ref, "N")
    table = [ref]
    for alt in (record.alts or ()):
        table.append(_normalize_allele(alt, ref))
    return table


def get_matrix_column(record, sample_names, resolve_IUPAC):
    """
    Transform one BCF record into a site-major nucleotide column.

    Unlike vcf2phylip's text parser, allele indices remain integers, so
    ALT indices >=10 are handled correctly rather than character-by-character.
    """
    table = _allele_table(record)
    column = []

    for name in sample_names:
        gt = record.samples[name].get("GT")
        if gt is None or len(gt) == 0:
            gt = (None,)

        try:
            if resolve_IUPAC:
                chosen = random.choice(gt)
                allele = "N" if chosen is None else table[chosen]
                key = "".join(sorted(set(allele)))
            else:
                nucs = []
                for idx in gt:
                    nucs.append("N" if idx is None else table[idx])
                key = "".join(sorted(set(nucs)))
        except (IndexError, TypeError):
            return "malformed"

        try:
            column.append(AMBIG[key])
        except KeyError:
            return "malformed"

    return "".join(column)


def get_matrix_column_bin(record, sample_names):
    """Return SNAPP-style 0/1/2/? states for diploid biallelic genotypes."""
    column = []
    for name in sample_names:
        gt = record.samples[name].get("GT")
        if gt is None:
            column.append("?")
        else:
            column.append(GEN_BIN_TUPLES.get(tuple(gt), "?"))
    return "".join(column)


def transpose_temp_to_sequence(temp_path, sample_index):
    """Read one site-major temporary matrix and extract one sample sequence."""
    chars = []
    with open(temp_path, "rt") as tmp_seq:
        for line in tmp_seq:
            line = line.rstrip("\n")
            chars.append(line[sample_index])
    return "".join(chars)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", action="store", dest="filename", required=True,
        help="Name of the input BCF file"
    )
    parser.add_argument(
        "--output-folder", action="store", dest="folder", default="./",
        help="Output folder name (current directory by default)"
    )
    parser.add_argument(
        "--output-prefix", action="store", dest="prefix",
        help="Prefix for output filenames (input BCF filename without extension by default)"
    )
    parser.add_argument(
        "-m", "--min-samples-locus", action="store",
        dest="min_samples_locus", type=int, default=4,
        help="Minimum samples required to be present at a locus (default=4)"
    )
    parser.add_argument(
        "-o", "--outgroup", action="store", dest="outgroup", default="",
        help="Name of outgroup. Sequence is written as first taxon."
    )
    parser.add_argument(
        "-p", "--phylip-disable", action="store_true", dest="phylipdisable",
        help="A PHYLIP matrix is written by default unless this flag is enabled"
    )
    parser.add_argument(
        "-f", "--fasta", action="store_true", dest="fasta",
        help="Write a FASTA matrix (disabled by default)"
    )
    parser.add_argument(
        "-n", "--nexus", action="store_true", dest="nexus",
        help="Write a NEXUS matrix (disabled by default)"
    )
    parser.add_argument(
        "-b", "--nexus-binary", action="store_true", dest="nexusbin",
        help="Write binary NEXUS for biallelic diploid SNPs (disabled by default)"
    )
    parser.add_argument(
        "-r", "--resolve-IUPAC", action="store_true", dest="resolve_IUPAC",
        help="Randomly resolve mixed genotypes instead of IUPAC ambiguities"
    )
    parser.add_argument(
        "-w", "--write-used-sites", action="store_true", dest="write_used",
        help="Save coordinates that passed filters"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Optional random seed for reproducible --resolve-IUPAC behavior"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s {version}".format(version=__version__)
    )
    args = parser.parse_args(argv)

    require_pysam()
    if args.seed is not None:
        random.seed(args.seed)

    outgroup = args.outgroup.split(",")[0].split(";")[0]

    input_path = Path(args.filename)
    if not input_path.exists():
        print("\nInput BCF file not found, please verify the provided path")
        return 1

    try:
        bcf = pysam.VariantFile(str(input_path), "r")
    except Exception as exc:
        print(f"\nCould not open BCF file: {exc}\n")
        return 1

    sample_names = extract_sample_names(bcf)
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in BCF; file may be corrupt or missing samples.\n")
        bcf.close()
        return 1

    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in BCF: {:d}".format(num_samples))

    args.min_samples_locus = min(num_samples, args.min_samples_locus)

    if not args.prefix:
        name = input_path.name
        if name.lower().endswith(".bcf"):
            name = name[:-4]
        args.prefix = name
    args.prefix += ".min" + str(args.min_samples_locus)

    outfolder = Path(args.folder)
    outfolder.mkdir(parents=True, exist_ok=True)
    outfile = str(outfolder / args.prefix)

    temporal = None
    temporalbin = None
    if args.fasta or args.nexus or not args.phylipdisable:
        temporal = open(outfile + ".tmp", "w")
    if args.nexusbin:
        temporalbin = open(outfile + ".bin.tmp", "w")

    used_sites = None
    if args.write_used:
        used_sites = open(outfile + ".used_sites.tsv", "w")
        used_sites.write("#CHROM\tPOS\tNUM_SAMPLES\n")

    snp_num = 0
    snp_accepted = 0
    snp_shallow = 0
    mnp_num = 0
    snp_biallelic = 0
    malformed_num = 0

    for record in bcf:
        snp_num += 1
        if snp_num % 500000 == 0:
            print("{:d} genotypes processed.".format(snp_num))

        num_samples_locus = num_genotypes(record, sample_names)
        if num_samples_locus < args.min_samples_locus:
            snp_shallow += 1
            continue

        if not is_snp(record):
            mnp_num += 1
            continue

        if args.fasta or args.nexus or not args.phylipdisable:
            site_tmp = get_matrix_column(record, sample_names, args.resolve_IUPAC)
            if site_tmp == "malformed":
                malformed_num += 1
                continue
            snp_accepted += 1
            temporal.write(site_tmp + "\n")
            if used_sites is not None:
                used_sites.write(
                    f"{record.contig}\t{record.pos}\t{num_samples_locus}\n"
                )

        if args.nexusbin:
            alts = record.alts or ()
            if len(alts) == 1 and len(_normalize_allele(alts[0], record.ref)) == 1:
                snp_biallelic += 1
                temporalbin.write(get_matrix_column_bin(record, sample_names) + "\n")

    bcf.close()

    print("Total of genotypes processed: {:d}".format(snp_num))
    print(
        "Genotypes excluded because they exceeded the amount of missing data allowed: "
        "{:d}".format(snp_shallow)
    )
    print(
        "Genotypes that passed missing data filter but were excluded for being MNPs: "
        "{:d}".format(mnp_num)
    )
    if malformed_num:
        print("Malformed SNP records excluded: {:d}".format(malformed_num))
    print("SNPs that passed the filters: {:d}".format(snp_accepted))
    if args.nexusbin:
        print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))

    if used_sites is not None:
        print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
        used_sites.close()

    if temporal is not None:
        temporal.close()
    if temporalbin is not None:
        temporalbin.close()

    output_phy = None
    output_fas = None
    output_nex = None
    output_nexbin = None

    if not args.phylipdisable:
        output_phy = open(outfile + ".phy", "w")
        output_phy.write("{:d} {:d}\n".format(len(sample_names), snp_accepted))

    if args.fasta:
        output_fas = open(outfile + ".fasta", "w")

    if args.nexus:
        output_nex = open(outfile + ".nexus", "w")
        output_nex.write(
            "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n"
            "\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n".format(
                len(sample_names), snp_accepted
            )
        )

    if args.nexusbin:
        output_nexbin = open(outfile + ".bin.nexus", "w")
        output_nexbin.write(
            "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n"
            "\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n".format(
                len(sample_names), snp_biallelic
            )
        )

    len_longest_name = max(len(name) for name in sample_names)
    idx_outgroup = sample_names.index(outgroup) if outgroup in sample_names else None

    def write_sample(index, label_prefix="Sample"):
        name = sample_names[index]
        padding = (len_longest_name + 3 - len(name)) * " "

        if temporal is not None:
            seqout = transpose_temp_to_sequence(outfile + ".tmp", index)
            if output_fas is not None:
                output_fas.write(">" + name + "\n" + seqout + "\n")
            if output_phy is not None:
                output_phy.write(name + padding + seqout + "\n")
            if output_nex is not None:
                output_nex.write(name + padding + seqout + "\n")

        if temporalbin is not None:
            seqbin = transpose_temp_to_sequence(outfile + ".bin.tmp", index)
            output_nexbin.write(name + padding + seqbin + "\n")

        if label_prefix == "Outgroup":
            print("Outgroup, '{}', added to the matrix(ces).".format(name))
        else:
            print(
                "Sample {:d} of {:d}, '{}', added to the matrix(ces).".format(
                    index + 1, len(sample_names), name
                )
            )

    if idx_outgroup is not None:
        write_sample(idx_outgroup, "Outgroup")

    for s in range(len(sample_names)):
        if s != idx_outgroup:
            write_sample(s)

    print()
    if output_phy is not None:
        print("PHYLIP matrix saved to: " + outfile + ".phy")
        output_phy.close()
    if output_fas is not None:
        print("FASTA matrix saved to: " + outfile + ".fasta")
        output_fas.close()
    if output_nex is not None:
        output_nex.write(";\nEND;\n")
        print("NEXUS matrix saved to: " + outfile + ".nexus")
        output_nex.close()
    if output_nexbin is not None:
        output_nexbin.write(";\nEND;\n")
        print("BINARY NEXUS matrix saved to: " + outfile + ".bin.nexus")
        output_nexbin.close()

    if temporal is not None:
        Path(outfile + ".tmp").unlink()
    if temporalbin is not None:
        Path(outfile + ".bin.tmp").unlink()

    print("\nDone!\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
