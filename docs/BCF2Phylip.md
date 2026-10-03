# BCF2Phylip

## Purpose

BCF2Phylip converts a multi-sample BCF directly into phylogenetic character
matrices without first materializing a VCF.

## Upstream reference

BCF2Phylip is derived from and benchmarked against VCF2Phylip:

https://github.com/edgardomortiz/vcf2phylip

https://doi.org/10.5281/zenodo.2540861

## Initial benchmark

The first benchmark should use a small jointly called haploid bacterial dataset,
preferably five *Mycobacterium tuberculosis* isolates.

Compare:

```text
BCF → BCF2Phylip → PHYLIP

BCF → bcftools view → VCF → VCF2Phylip → PHYLIP
```

The primary correctness target is equality of:

- retained sites;
- sample order;
- alignment length;
- character state for every sample at every retained site.

## Later performance benchmark

Compare:

1. direct BCF → BCF2Phylip;
2. BCF → plain VCF → VCF2Phylip;
3. BCF → VCF.gz → VCF2Phylip.

Measure:

- wall time;
- peak RAM;
- temporary disk use;
- total input/output volume;
- final matrix identity.
