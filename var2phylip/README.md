# Future importable Python package

This directory is reserved for reusable VAR2Phylip modules.

The initial prototype scripts should first be placed under:

```text
scripts/BCF2Phylip.py
scripts/PAF2Phylip.py
```

After validation, shared functionality can be refactored here.

Possible future layout:

```text
var2phylip/
├── __init__.py
├── bcf.py
├── paf.py
├── matrix.py
├── iupac.py
└── validation.py
```
