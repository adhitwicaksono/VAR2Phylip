# Credits and upstream software

## VCF2Phylip

BCF2Phylip builds upon the design and genotype-to-matrix behavior of
**VCF2Phylip** by Edgardo M. Ortiz, with credits to Juan D. Palacio-Mejía.

Original repository:

https://github.com/edgardomortiz/vcf2phylip

Zenodo software record:

https://doi.org/10.5281/zenodo.2540861

VCF2Phylip is licensed under GPL-3.0.

The original VCF2Phylip program is not redistributed inside VAR2Phylip.
BCF2Phylip source files should retain explicit upstream attribution and
modification notices.

## PAF2Phylip

PAF2Phylip is an independent implementation developed for VAR2Phylip.

Its design relies on the PAF format and `cs` tag produced by minimap2, but
minimap2 source code is not incorporated into PAF2Phylip.

Users should cite minimap2 when appropriate for analyses that use minimap2-
generated alignments.
