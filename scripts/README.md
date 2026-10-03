# Scripts

Add the current prototypes here:

```text
scripts/BCF2Phylip.py
scripts/PAF2Phylip.py
```

Recommended next step after the first public benchmark:

```text
scripts/
    BCF2Phylip.py       lightweight CLI wrapper
    PAF2Phylip.py       lightweight CLI wrapper

var2phylip/
    bcf.py              reusable BCF logic
    paf.py              reusable PAF logic
    matrix.py           shared output functions
    validation.py       benchmark/comparison utilities
```

Do not move code into the package module until the prototype behavior is frozen
and regression-tested.
