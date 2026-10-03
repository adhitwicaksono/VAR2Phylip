# VAR2Phylip

**VAR2Phylip** is a small bioinformatics toolkit for converting compact genomic
variant and alignment representations directly into phylogenetic character
matrices.

The project currently contains two planned tools:

- **BCF2Phylip** — direct BCF → PHYLIP / FASTA / NEXUS conversion without
  materializing an intermediate VCF.
- **PAF2Phylip** — reference-anchored SNP matrices from minimap2 PAF
  alignments containing `cs` tags.

The project is initially focused on **closely related haploid prokaryotic
genomes**, where direct benchmarking and manual validation are tractable.

---

## Relationship to VCF2Phylip

BCF2Phylip is derived from and intended to preserve the relevant
genotype-to-matrix behavior of **VCF2Phylip** by Edgardo M. Ortiz, with credits
to Juan D. Palacio-Mejía.

Original project:

https://github.com/edgardomortiz/vcf2phylip

Zenodo software record:

https://doi.org/10.5281/zenodo.2540861

VCF2Phylip itself is **not redistributed** in this repository.

Because BCF2Phylip is a derivative work of GPL-3.0-licensed software,
VAR2Phylip is distributed under the **GNU General Public License v3.0**.

---

## Project concept

VAR2Phylip asks a simple question:

> Can machine-efficient genomic representations be converted directly into
> phylogenetic matrices without unnecessary intermediate formats?

The two current routes are:

```text
BCF
 ↓
structured genotypes
 ↓
SNP characters
 ↓
PHYLIP / FASTA / NEXUS
```

and:

```text
PAF + reference
 ↓
reference-coordinate sequence states
 ↓
SNP characters
 ↓
PHYLIP
```

---

## Repository layout

```text
VAR2Phylip/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
├── CREDITS.md
├── pyproject.toml
├── .gitignore
│
├── scripts/
│   ├── BCF2Phylip.py        # add prototype here
│   └── PAF2Phylip.py        # add prototype here
│
├── var2phylip/
│   └── README.md             # future importable package layout
│
├── docs/
│   ├── BCF2Phylip.md
│   ├── PAF2Phylip.md
│   └── VCF2Phylip_reference.md
│
├── benchmark/
│   ├── README.md
│   ├── MTB_BCF_benchmark.md
│   └── MTB_PAF_benchmark.md
│
└── tests/
    └── README.md
```

---

## Add the prototype scripts

This initial repository skeleton intentionally does **not** include the Python
prototype files.

Add them later as:

```text
scripts/BCF2Phylip.py
scripts/PAF2Phylip.py
```

Once the prototypes stabilize, their reusable functions can be migrated into:

```text
var2phylip/
```

while the files under `scripts/` can remain lightweight command-line entry
points.

---

## BCF2Phylip

### Goal

Read multi-sample BCF directly through HTSlib/pysam and produce the same
phylogenetic character matrix expected from an equivalent VCF processed by
VCF2Phylip.

The initial benchmark contract is:

```text
same multi-sample BCF
       │
       ├── BCF2Phylip → direct.phy
       │
       └── bcftools view → VCF → VCF2Phylip → reference.phy

direct.phy == reference.phy
```

The first validation target is a small jointly called **Mycobacterium
tuberculosis** dataset.

### Initial compatibility targets

- SNP-only filtering
- multiallelic SNP support
- haploid and arbitrary-ploidy genotype handling
- IUPAC consensus
- missing-data thresholds
- `<NON_REF>` handling
- `*` deletion handling
- sample order
- outgroup-first output
- used-site tracking

The prototype intentionally preserves relevant VCF2Phylip behavior while using
structured BCF genotype access instead of text-field parsing.

---

## PAF2Phylip

### Goal

Build a conservative reference-anchored SNP matrix directly from multiple
assembly-to-reference PAF alignments generated with minimap2 `--cs`.

Example conceptual workflow:

```text
reference.fa
   │
   ├── sample01.fa → minimap2 --cs → sample01.paf
   ├── sample02.fa → minimap2 --cs → sample02.paf
   ├── sample03.fa → minimap2 --cs → sample03.paf
   └── sample04.fa → minimap2 --cs → sample04.paf
                                   ↓
                              PAF2Phylip
                                   ↓
                         reference-anchored SNPs
                                   ↓
                                PHYLIP
```

Version 0.1 is intentionally limited to closely related haploid genomes and
does not attempt to build a general structural-variant-aware whole-genome
multiple alignment.

---

## Benchmark philosophy

Correctness comes before optimization.

The planned order is:

1. synthetic edge-case fixtures;
2. small real MTB dataset;
3. exact site and matrix comparison;
4. larger bacterial cohort;
5. runtime, RAM, and temporary-disk benchmarking;
6. only then, optimization of matrix transposition and I/O.

For BCF2Phylip, VCF2Phylip is the primary behavioral reference.

For PAF2Phylip, the primary validation route is direct comparison with minimap2
`cs` information and paftools-derived calls under controlled callability
semantics.

---

## Status

**Prototype stage / pre-release.**

Nothing in this repository should yet be treated as a validated production
bioinformatics workflow.

Suggested versioning:

```text
0.1.0-alpha   first public prototype
0.2.0-alpha   first real-data benchmark
0.3.0-beta    stabilized CLI and regression suite
1.0.0         validated initial release
```

---

## Citation

Until VAR2Phylip receives its own DOI or publication, cite the repository
version/commit used in your analysis.

For BCF2Phylip, also cite the upstream VCF2Phylip software:

> Ortiz EM. VCF2Phylip. Zenodo. https://doi.org/10.5281/zenodo.2540861

See `CITATION.cff` and `CREDITS.md`.

---

## License

GNU General Public License v3.0.

See `LICENSE`.
