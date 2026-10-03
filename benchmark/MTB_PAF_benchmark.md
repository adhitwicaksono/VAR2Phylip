# MTB PAF benchmark

## Goal

Validate PAF2Phylip on a small collection of closely related haploid bacterial
assemblies.

## Proposed workflow

```text
reference assembly
       │
       ├── isolate01 assembly → minimap2 --cs → isolate01.paf
       ├── isolate02 assembly → minimap2 --cs → isolate02.paf
       ├── isolate03 assembly → minimap2 --cs → isolate03.paf
       ├── isolate04 assembly → minimap2 --cs → isolate04.paf
       └── isolate05 assembly → minimap2 --cs → isolate05.paf
```

## Validation

For every retained PAF2Phylip SNP:

- confirm reference coordinate;
- confirm REF base;
- confirm sample-specific ALT state;
- confirm callability;
- compare against controlled paftools-derived calls.

The benchmark should explicitly document any difference in callability policy.
