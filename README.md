# Topology-Spectral-Sequence

Compute and inspect algebraic spectral sequences over SymPy domains, then export chart data to HTML.

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
p3 = ss.add_page({t: a, 2*a: 0})
p4 = ss.add_page()

module = p4[3, 2]          # module at bidegree (3, 2)
info = module.get_structural_information()
if info is not None:
    gens, torsion = info
```

### Important API notes

- `ss.add_page(known_diff)` expects a dict of SymPy expressions in your generators.
- `p = ss.add_page(...)` returns a `Page`; index modules with `p[x, y]`.
- `module.get_structural_information()` returns `(generators, torsion)` in absolute coordinates.
- `module.get_diff_span()` and `page.d.get_diff_span(bidegree)` compute differential images. If data is insufficient, the program may request missing differential values interactively.

## Base-domain assumptions (read this first)

This code currently assumes a Euclidean-style computational domain for Smith normal form routines.

- Internally, SNF uses operations like `gcdex`, `rem`, and `exquo`.
- `SNFMatrix` also verifies `domain.is_PID`.
- In practice, use Euclidean domains (for example `ZZ`, finite fields like `GF(p)`, and other fields with the required exact operations).
- A non-Euclidean PID is not a supported target here; decomposition can fail or hit the explicit non-convergence guard.

If you provide a custom domain, ensure those exact-division and gcd/remainder operations are implemented and consistent.
