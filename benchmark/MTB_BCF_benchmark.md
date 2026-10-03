# MTB BCF benchmark

## Goal

Show that direct BCF parsing produces the same phylogenetic matrix as the
equivalent VCF processed by VCF2Phylip.

## Proposed workflow

```text
5 MTB FASTQ datasets
        ↓
same mapping and calling pipeline
        ↓
joint multi-sample BCF
       ↙        ↘
BCF2Phylip     bcftools view
     ↓              ↓
direct.phy          VCF
                     ↓
               VCF2Phylip
                     ↓
               reference.phy
```

## First comparison

Run strict presence:

```text
-m 5
```

Then repeat:

```text
-m 4
```

to exercise missing-data behavior.

## Record

- source accessions;
- reference genome;
- aligner and version;
- variant caller and version;
- exact command lines;
- BCF size;
- VCF size;
- VCF.gz size;
- retained sites;
- output hashes;
- wall time;
- peak RAM;
- temporary disk use.
