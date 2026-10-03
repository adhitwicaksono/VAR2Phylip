#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAF2Phylip v0.1

Build a reference-anchored SNP PHYLIP matrix directly from multiple
assembly-to-reference minimap2 PAF files carrying the cs:Z tag.

Designed first for closely related haploid prokaryotic assemblies.

Core idea
---------
reference FASTA + one PAF per sample
              |
              v
  shared reference-coordinate SNP sites
              |
              v
       relaxed PHYLIP matrix

This prototype intentionally does NOT try to build a full gapped multiple
whole-genome alignment. Insertions relative to the reference are ignored in
v0.1; deletions make the affected reference coordinates uncallable for that
sample.

PAF2Phylip is an independent implementation. No minimap2/paftools source code
is incorporated.
"""

from __future__ import annotations

import argparse
import bisect
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

__author__ = "Adhityo Wicaksono"
__version__ = "0.1.0"
__date__ = "2026-10-04"

CS_TOKEN_RE = re.compile(
    r":[0-9]+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|[+\-][A-Za-z]+|~[A-Za-z]{2}[0-9]+[A-Za-z]{2}"
)

DNA = {"A", "C", "G", "T"}


@dataclass
class IntervalIndex:
    """Count half-open interval coverage at a coordinate."""

    starts: Dict[str, List[int]] = field(default_factory=lambda: defaultdict(list))
    ends: Dict[str, List[int]] = field(default_factory=lambda: defaultdict(list))

    def add(self, contig: str, start: int, end: int) -> None:
        if end <= start:
            return
        self.starts[contig].append(start)
        self.ends[contig].append(end)

    def finalize(self) -> None:
        for contig in set(self.starts) | set(self.ends):
            self.starts[contig].sort()
            self.ends[contig].sort()

    def count(self, contig: str, pos0: int) -> int:
        starts = self.starts.get(contig, ())
        ends = self.ends.get(contig, ())
        # active half-open intervals [start, end):
        # starts <= pos minus ends <= pos
        return bisect.bisect_right(starts, pos0) - bisect.bisect_right(ends, pos0)


@dataclass
class SampleModel:
    sample: str
    paf_path: Path
    coverage: IntervalIndex = field(default_factory=IntervalIndex)
    deletions: IntervalIndex = field(default_factory=IntervalIndex)
    substitutions: Dict[Tuple[str, int], Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    alignments_seen: int = 0
    alignments_used: int = 0
    skipped_secondary: int = 0
    skipped_low_mapq: int = 0
    skipped_short: int = 0
    substitution_events: int = 0
    deletion_bases: int = 0
    insertion_bases: int = 0

    def finalize(self) -> None:
        self.coverage.finalize()
        self.deletions.finalize()


def parse_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name: Optional[str] = None

    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name:
                    raise ValueError("FASTA contains an empty sequence name")
                if name in seqs:
                    raise ValueError(f"Duplicate FASTA sequence name: {name}")
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence found before first header")
                seqs[name].append(line.upper())

    if not seqs:
        raise ValueError(f"No sequences found in reference FASTA: {path}")

    joined = {name: "".join(parts) for name, parts in seqs.items()}
    for name, seq in joined.items():
        if not seq:
            raise ValueError(f"Reference sequence is empty: {name}")
    return joined


def parse_sample_sheet(path: Path) -> List[Tuple[str, Path]]:
    rows: List[Tuple[str, Path]] = []
    seen: Set[str] = set()
    base = path.parent

    with path.open("rt", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_no}: expected exactly 2 tab-separated columns: "
                    "SAMPLE<TAB>PAF"
                )
            sample, paf_text = fields
            if any(ch.isspace() for ch in sample):
                raise ValueError(
                    f"{path}:{line_no}: sample names may not contain whitespace: {sample!r}"
                )
            if sample in seen:
                raise ValueError(f"{path}:{line_no}: duplicate sample name: {sample}")
            seen.add(sample)

            paf = Path(paf_text)
            if not paf.is_absolute():
                paf = (base / paf).resolve()
            if not paf.is_file():
                raise FileNotFoundError(f"{path}:{line_no}: PAF not found: {paf}")

            rows.append((sample, paf))

    if not rows:
        raise ValueError(f"No samples found in sample sheet: {path}")
    return rows


def parse_tags(fields: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    tags: Dict[str, Tuple[str, str]] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tag, typ, value = parts
            tags[tag] = (typ, value)
    return tags


def parse_cs(
    cs: str,
    contig: str,
    target_start: int,
    reference: Dict[str, str],
    model: SampleModel,
    source: str,
) -> None:
    tokens = CS_TOKEN_RE.findall(cs)
    if "".join(tokens) != cs:
        raise ValueError(f"{source}: unsupported or malformed cs tag: {cs}")

    if contig not in reference:
        raise ValueError(
            f"{source}: target contig {contig!r} is absent from reference FASTA"
        )

    refseq = reference[contig]
    rpos = target_start

    for token in tokens:
        op = token[0]

        if op == ":":
            rpos += int(token[1:])

        elif op == "=":
            seq = token[1:].upper()
            expected = refseq[rpos:rpos + len(seq)]
            if expected != seq:
                raise ValueError(
                    f"{source}: cs long-match disagrees with reference at "
                    f"{contig}:{rpos + 1}. cs={seq!r}, reference={expected!r}"
                )
            rpos += len(seq)

        elif op == "*":
            ref_base = token[1].upper()
            qry_base = token[2].upper()
            if rpos >= len(refseq):
                raise ValueError(f"{source}: cs substitution exceeds reference length")
            expected = refseq[rpos]
            if expected != ref_base:
                raise ValueError(
                    f"{source}: cs substitution disagrees with reference at "
                    f"{contig}:{rpos + 1}. cs REF={ref_base}, FASTA REF={expected}"
                )
            model.substitutions[(contig, rpos)].add(qry_base)
            model.substitution_events += 1
            rpos += 1

        elif op == "+":
            # Query insertion relative to target. It has no reference coordinate
            # and is intentionally ignored in SNP-only v0.1.
            model.insertion_bases += len(token) - 1

        elif op == "-":
            deleted = token[1:].upper()
            expected = refseq[rpos:rpos + len(deleted)]
            if expected != deleted:
                raise ValueError(
                    f"{source}: cs deletion disagrees with reference at "
                    f"{contig}:{rpos + 1}. cs={deleted!r}, reference={expected!r}"
                )
            model.deletions.add(contig, rpos, rpos + len(deleted))
            model.deletion_bases += len(deleted)
            rpos += len(deleted)

        elif op == "~":
            raise ValueError(
                f"{source}: splice/intron cs operation {token!r} is unsupported in "
                "PAF2Phylip v0.1. Use genomic assembly-to-reference PAFs."
            )

        else:
            raise AssertionError(f"Unhandled cs operation: {token}")

    return


def load_paf(
    sample: str,
    paf_path: Path,
    reference: Dict[str, str],
    min_mapq: int,
    min_aln_len: int,
    include_secondary: bool,
) -> SampleModel:
    model = SampleModel(sample=sample, paf_path=paf_path)

    with paf_path.open("rt", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"{paf_path}:{line_no}: expected >=12 PAF columns, found {len(fields)}"
                )

            model.alignments_seen += 1

            try:
                qname = fields[0]
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                strand = fields[4]
                tname = fields[5]
                tlen = int(fields[6])
                tstart = int(fields[7])
                tend = int(fields[8])
                nmatch = int(fields[9])
                block_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError as exc:
                raise ValueError(f"{paf_path}:{line_no}: malformed numeric PAF field") from exc

            if strand not in {"+", "-"}:
                raise ValueError(f"{paf_path}:{line_no}: invalid strand: {strand}")

            if tname not in reference:
                raise ValueError(
                    f"{paf_path}:{line_no}: target {tname!r} not found in reference FASTA"
                )
            if tlen != len(reference[tname]):
                raise ValueError(
                    f"{paf_path}:{line_no}: PAF target length {tlen} differs from "
                    f"reference FASTA length {len(reference[tname])} for {tname}"
                )
            if not (0 <= tstart <= tend <= tlen):
                raise ValueError(f"{paf_path}:{line_no}: invalid target coordinates")

            tags = parse_tags(fields[12:])
            tp = tags.get("tp")
            if not include_secondary and tp is not None and tp[1] == "S":
                model.skipped_secondary += 1
                continue

            if mapq < min_mapq:
                model.skipped_low_mapq += 1
                continue

            if (tend - tstart) < min_aln_len:
                model.skipped_short += 1
                continue

            cs_tag = tags.get("cs")
            if cs_tag is None or cs_tag[0] != "Z":
                raise ValueError(
                    f"{paf_path}:{line_no}: missing cs:Z tag. Generate PAF with minimap2 --cs."
                )

            model.coverage.add(tname, tstart, tend)
            parse_cs(
                cs=cs_tag[1],
                contig=tname,
                target_start=tstart,
                reference=reference,
                model=model,
                source=f"{paf_path}:{line_no}",
            )
            model.alignments_used += 1

    model.finalize()

    if model.alignments_used == 0:
        raise ValueError(f"No usable alignments remained for sample {sample}: {paf_path}")

    return model


def sample_state(
    model: SampleModel,
    contig: str,
    pos0: int,
    ref_base: str,
) -> Optional[str]:
    """
    Return a callable nucleotide state, otherwise None.

    A site is callable only if exactly one accepted PAF alignment covers the
    reference coordinate and the coordinate is not deleted in that sample.
    This deliberately rejects overlap/duplication ambiguity.
    """
    if model.coverage.count(contig, pos0) != 1:
        return None
    if model.deletions.count(contig, pos0) > 0:
        return None

    alts = model.substitutions.get((contig, pos0))
    if not alts:
        return ref_base
    if len(alts) != 1:
        return None

    alt = next(iter(alts)).upper()
    if alt not in DNA:
        return None
    return alt


def write_phylip(
    path: Path,
    names: Sequence[str],
    sequences: Dict[str, str],
    nchar: int,
) -> None:
    longest = max(len(name) for name in names)
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(names)} {nchar}\n")
        for name in names:
            seq = sequences[name]
            if len(seq) != nchar:
                raise RuntimeError(
                    f"Internal sequence-length mismatch for {name}: {len(seq)} != {nchar}"
                )
            handle.write(f"{name}{' ' * (longest + 3 - len(name))}{seq}\n")


def write_fasta(path: Path, names: Sequence[str], sequences: Dict[str, str]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for name in names:
            handle.write(f">{name}\n{sequences[name]}\n")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Convert multiple assembly-to-reference minimap2 PAFs with cs tags "
            "directly into a reference-anchored SNP PHYLIP matrix."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-r", "--reference", required=True, help="Reference FASTA used as minimap2 target")
    ap.add_argument(
        "-s", "--samples", required=True,
        help="TSV with two columns: SAMPLE<TAB>PAF"
    )
    ap.add_argument(
        "-o", "--output-prefix", default="PAF2Phylip",
        help="Output path prefix"
    )
    ap.add_argument(
        "-m", "--min-samples-locus", type=int, default=None,
        help="Minimum callable query samples at a SNP; default = all samples"
    )
    ap.add_argument(
        "--min-mapq", type=int, default=0,
        help="Minimum PAF mapping quality"
    )
    ap.add_argument(
        "--min-aln-len", type=int, default=0,
        help="Minimum accepted target-alignment span"
    )
    ap.add_argument(
        "--include-secondary", action="store_true",
        help="Include alignments tagged tp:A:S (not recommended for core-SNP inference)"
    )
    ap.add_argument(
        "--include-reference", action="store_true",
        help="Include the reference as an additional taxon in the output matrix"
    )
    ap.add_argument(
        "--reference-name", default="REFERENCE",
        help="Taxon label used with --include-reference"
    )
    ap.add_argument(
        "--fasta", action="store_true",
        help="Also write FASTA alignment"
    )
    ap.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    ref_path = Path(args.reference).resolve()
    sheet_path = Path(args.samples).resolve()
    prefix = Path(args.output_prefix)

    if not ref_path.is_file():
        raise SystemExit(f"Reference FASTA not found: {ref_path}")
    if not sheet_path.is_file():
        raise SystemExit(f"Sample sheet not found: {sheet_path}")

    try:
        reference = parse_fasta(ref_path)
        sample_rows = parse_sample_sheet(sheet_path)

        models: List[SampleModel] = []
        for sample, paf in sample_rows:
            models.append(
                load_paf(
                    sample=sample,
                    paf_path=paf,
                    reference=reference,
                    min_mapq=args.min_mapq,
                    min_aln_len=args.min_aln_len,
                    include_secondary=args.include_secondary,
                )
            )
    except (ValueError, FileNotFoundError) as exc:
        raise SystemExit(str(exc)) from exc

    query_names = [model.sample for model in models]
    n_queries = len(query_names)

    min_samples = args.min_samples_locus
    if min_samples is None:
        min_samples = n_queries
    if min_samples < 1:
        raise SystemExit("--min-samples-locus must be >=1")
    if min_samples > n_queries:
        min_samples = n_queries

    # Candidate positions are the union of substitutions observed in accepted PAFs.
    candidates: Set[Tuple[str, int]] = set()
    for model in models:
        candidates.update(model.substitutions.keys())

    contig_order = {name: i for i, name in enumerate(reference)}
    ordered_candidates = sorted(
        candidates,
        key=lambda key: (contig_order[key[0]], key[1]),
    )

    seq_chars: Dict[str, List[str]] = {name: [] for name in query_names}
    if args.include_reference:
        if args.reference_name in seq_chars:
            raise SystemExit(
                f"--reference-name conflicts with query sample name: {args.reference_name}"
            )
        seq_chars[args.reference_name] = []

    used_rows: List[Tuple[str, int, str, int, str]] = []
    excluded_low_presence = 0
    excluded_monomorphic = 0

    for contig, pos0 in ordered_candidates:
        ref_base = reference[contig][pos0].upper()
        states: List[Optional[str]] = [
            sample_state(model, contig, pos0, ref_base) for model in models
        ]
        called = sum(state is not None for state in states)

        if called < min_samples:
            excluded_low_presence += 1
            continue

        observed = {state for state in states if state is not None}
        if len(observed) < 2:
            excluded_monomorphic += 1
            continue

        for name, state in zip(query_names, states):
            seq_chars[name].append(state if state is not None else "N")

        if args.include_reference:
            seq_chars[args.reference_name].append(ref_base)

        used_rows.append(
            (
                contig,
                pos0 + 1,
                ref_base,
                called,
                ",".join(sorted(observed)),
            )
        )

    names = list(query_names)
    if args.include_reference:
        names.insert(0, args.reference_name)

    sequences = {name: "".join(chars) for name, chars in seq_chars.items()}
    nchar = len(used_rows)

    phy_path = prefix.with_suffix(".phy")
    used_path = prefix.with_suffix(".used_sites.tsv")
    summary_path = prefix.with_suffix(".summary.tsv")
    fasta_path = prefix.with_suffix(".fasta")

    phy_path.parent.mkdir(parents=True, exist_ok=True)
    write_phylip(phy_path, names, sequences, nchar)

    with used_path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("#CHROM\tPOS\tREF\tCALLED_SAMPLES\tOBSERVED_STATES\n")
        for row in used_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("sample\tpaf\talignments_seen\talignments_used\t"
                     "skipped_secondary\tskipped_low_mapq\tskipped_short\t"
                     "substitution_events\tdeletion_bases\tinsertion_bases\n")
        for model in models:
            handle.write(
                f"{model.sample}\t{model.paf_path}\t{model.alignments_seen}\t"
                f"{model.alignments_used}\t{model.skipped_secondary}\t"
                f"{model.skipped_low_mapq}\t{model.skipped_short}\t"
                f"{model.substitution_events}\t{model.deletion_bases}\t"
                f"{model.insertion_bases}\n"
            )

    if args.fasta:
        write_fasta(fasta_path, names, sequences)

    print(f"Samples: {n_queries}")
    print(f"Candidate substitution sites: {len(ordered_candidates)}")
    print(f"Minimum callable samples per site: {min_samples}")
    print(f"Excluded for low presence: {excluded_low_presence}")
    print(f"Excluded as monomorphic after callability filtering: {excluded_monomorphic}")
    print(f"Accepted SNP sites: {nchar}")
    print(f"PHYLIP: {phy_path}")
    print(f"Used sites: {used_path}")
    print(f"Summary: {summary_path}")
    if args.fasta:
        print(f"FASTA: {fasta_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
