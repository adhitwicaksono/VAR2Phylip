# PAF2Phylip

## Purpose

PAF2Phylip builds conservative reference-anchored SNP matrices directly from
multiple assembly-to-reference PAF alignments carrying minimap2 `cs` tags.

## Initial scope

- haploid genomes;
- closely related prokaryotic assemblies;
- one common reference;
- one PAF per sample;
- SNP-only reference-coordinate characters;
- no general structural-variant-aware multiple alignment.

## Initial benchmark

Generate:

```bash
minimap2 -cx asm5 --cs reference.fa sample.fa > sample.paf
```

for each isolate.

Then compare accepted substitutions and callable positions against controlled
paftools-derived calls.

The important question is not merely whether substitutions match, but whether
callability semantics match.

## Later extensions

Possible later work:

- invariant core-genome alignment mode;
- explicit structural-variation exclusion reports;
- alternative callability thresholds;
- blockwise output for larger cohorts.

Insertions relative to the reference should not be incorporated into a
multi-genome alignment until a defensible cross-sample homology policy is
implemented.
