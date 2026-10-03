# Contributing

VAR2Phylip is currently an early-stage research software project.

Contributions are welcome, especially for:

- reproducible test datasets;
- genotype edge cases;
- BCF parsing and HTSlib behavior;
- PAF `cs` parsing;
- bacterial benchmark datasets;
- PHYLIP / FASTA / NEXUS compatibility;
- memory and disk-efficiency profiling;
- regression testing.

## Development principles

1. **Correctness before speed.**
2. **Keep benchmark semantics explicit.**
3. **Do not silently reinterpret missing data.**
4. **Preserve exact sample order unless the user explicitly requests otherwise.**
5. **Document any intentional behavior difference from VCF2Phylip.**
6. **Use small synthetic fixtures before large real datasets.**
7. **Avoid organism-specific assumptions unless the mode explicitly requires them.**

## Pull requests

Please include:

- a concise description of the change;
- the biological or technical rationale;
- a minimal test or fixture when possible;
- expected output;
- any compatibility impact.

For BCF2Phylip changes, note whether the behavior remains compatible with
VCF2Phylip and whether the difference is intentional.

## Style

Prefer readable Python over clever Python.

Public functions should have short docstrings. New command-line options should
be documented in both `--help` and the relevant file under `docs/`.
