# VCF2Phylip as the BCF2Phylip reference implementation

BCF2Phylip uses VCF2Phylip as its primary behavioral reference.

Original software:

**VCF2Phylip**  
Edgardo M. Ortiz  
with credits to Juan D. Palacio-Mejía

Repository:

https://github.com/edgardomortiz/vcf2phylip

Zenodo:

https://doi.org/10.5281/zenodo.2540861

## Compatibility principle

BCF2Phylip should preserve intended biological semantics where appropriate:

- SNP filtering;
- multiallelic SNP support;
- genotype-to-base conversion;
- IUPAC consensus;
- missing-data thresholds;
- sample ordering;
- used-site reporting.

BCF2Phylip does not need to preserve accidental limitations caused by
line-oriented VCF text parsing when structured BCF/HTSlib access provides a
more reliable representation.

Any intentional difference must be documented and tested.
