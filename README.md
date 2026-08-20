# Topology-Spectral-Sequence

Compute and inspect algebraic spectral sequences over exact SymPy domains, then export chart data to HTML.

## Run the program

From the repository root, run:

```bash
python src/main.py <min_x> <max_x> <min_y> <max_y> <out.html>
```

Example:

```bash
python src/main.py 0 8 0 8 output_current_main.html
```

This builds the sample spectral sequence in `src/main.py`, scans the rectangle of bidegrees you provide, and writes an HTML chart file.

## API for custom spectral sequences

The core API is in `src/spectral_sequence.py`.

```python
from sympy import ZZ
from sympy.abc import a, t
from src.spectral_sequence import SpectralSequence

ss = SpectralSequence(
    ZZ,                     # base domain
    [a, t],                 # generators
    [[3, 0], [0, 2]],       # generator bidegrees (2 x n)
    [[1, 0], [-1, 1]],      # differential-bidegree coefficient matrix
)

ss.kill(a**2)              # add relations

# Add pages (E1, E2, E3, ...). known_diff maps src -> d(src) on that page.
p1 = ss.add_page({a: 0, t: 0})
p2 = ss.add_page({a: 0, t: 0})
p3 = ss.add_page({t: a, a: 0})
p4 = ss.add_page()

module = p4[3, 2]          # module at bidegree (3, 2)
info = module.get_structural_information()
if info is not None:
    gens, torsion = info
```

In this sample, the supplied value `d3(a) = 0` is mathematically redundant because `d3(t) = a` and `d3^2 = 0`; the current implementation does not yet infer `d^2 = 0` automatically.

### Important API notes

- `ss.add_page(known_diff)` expects a dict of SymPy expressions in your generators.
- `p = ss.add_page(...)` returns a `Page`; index modules with `p[x, y]`.
- `module.get_structural_information()` returns `(generators, torsion)` in absolute coordinates.
- `module.get_diff_span()` and `page.d.get_diff_span(bidegree)` compute differential images. If data is insufficient, the program may request missing differential values interactively.

## Base-domain assumptions

This repository uses a **custom Smith-normal-form implementation** in `src/snf.py`, built on SymPy's exact `DomainMatrix`/domain arithmetic. It does not simply delegate the spectral-sequence computations to SymPy's high-level Smith-normal-form routine.

- The SNF implementation uses exact domain operations such as `gcdex`, `rem`, and `exquo`.
- `SNFMatrix` verifies `domain.is_PID`.
- In practice, the implementation is aimed at Euclidean-style computational PID domains such as `ZZ`, finite fields like `GF(p)`, and other exact fields/domains that provide the required operations.
- A non-Euclidean PID is not currently a supported computational target; decomposition can fail or hit the explicit non-convergence guard.

If you provide a custom domain, ensure those exact-division, gcd, and remainder operations are implemented consistently.
