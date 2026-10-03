# Benchmarks

Benchmarks serve two separate goals:

1. **correctness**
2. **performance**

Correctness benchmarks must be completed before performance claims are made.

## Planned benchmark tiers

### Tier 1 — synthetic fixtures

Exercise:

- haploid REF/ALT;
- missing genotypes;
- multiallelic SNPs;
- `<NON_REF>`;
- `*` deletion allele;
- arbitrary ploidy;
- PAF substitutions;
- PAF deletions;
- ambiguous/overlapping PAF coverage.

### Tier 2 — five-isolate MTB dataset

Use the same biological material to generate:

- a multi-sample BCF;
- a VCF derived from that BCF;
- assembly-to-reference PAFs when assemblies are available.

### Tier 3 — larger prokaryotic cohort

Measure scaling across sample count and SNP count.

### Tier 4 — large eukaryotic stress test

Only after the prokaryotic benchmark suite is stable.
