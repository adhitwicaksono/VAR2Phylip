# Tests

This directory is reserved for regression tests and small curated fixtures.

Recommended future structure:

```text
tests/
├── fixtures/
│   ├── bcf/
│   └── paf/
├── expected/
├── test_bcf_semantics.py
├── test_paf_cs.py
└── test_matrix_outputs.py
```

Large real datasets should not be committed to Git.

Small synthetic BCF/VCF/PAF/FASTA fixtures may be committed under
`tests/fixtures/`.
