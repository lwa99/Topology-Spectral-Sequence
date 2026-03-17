"""Examples for a cohomological spectral sequence project.

These examples are written to match the *screenshot-style* page bookkeeping used in the
conversation:

- `p_1` is an extra startup page,
- `p_2` is the mathematically natural input page for Serre-type examples,
- `p_r` is the page on which the differential `d_r` is specified,
- `p_{r+1}` is the next page obtained from that differential.

So this file is aligned with the user's current code style rather than the clean textbook
convention of starting directly at `E_2`.

Important notes:
1. These examples assume a cohomological grading convention
       d_r : E_r^{p,q} -> E_r^{p+r, q-r+1}.
2. The current project discussed in the chat does *not* automatically extend differentials
   by the Leibniz rule. Therefore, examples whose mathematical meaning depends on Leibniz
   should be treated as conceptual test cases unless the relevant values are supplied by hand.
3. The import path for `SpectralSequence` is project-specific. Replace the placeholder import
   below with the correct one in your codebase.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List

from sympy import symbols
from sympy.polys.domains import GF

# Example: from my_project.sseq import SpectralSequence
from src.spectral_sequence import SpectralSequence  # type: ignore


@dataclass(frozen=True)
class ExampleRecord:
    """Container for one example and its builder."""

    name: str
    description: str
    build: Callable[[], Dict[str, object]]


def build_example_1_trivial_product() -> Dict[str, object]:
    """Build the spectral sequence for S^1 -> S^2 x S^1 -> S^2 over F_2."""
    a, u = symbols("a u")
    ss = SpectralSequence(
        GF(2),
        [a, u],
        [[2, 0], [0, 1]],
        [[1, 0], [-1, 1]],
    )
    ss.kill(a**2)
    ss.kill(u**2)

    p_1 = ss.add_page({a: 0, u: 0})
    p_2 = ss.add_page({a: 0, u: 0})
    p_3 = ss.add_page({a: 0, u: 0})

    return {
        "name": "Trivial product fibration S^1 -> S^2 x S^1 -> S^2",
        "field": "GF(2)",
        "generators": {"a": a, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3},
        "ss": ss,
        "first_nonzero_differential": None,
        "expected_stable_page": "p_2",
    }


def build_example_2_hopf_fibration() -> Dict[str, object]:
    """Build the spectral sequence for the Hopf fibration S^1 -> S^3 -> S^2."""
    a, u = symbols("a u")
    ss = SpectralSequence(
        GF(2),
        [a, u],
        [[2, 0], [0, 1]],
        [[1, 0], [-1, 1]],
    )
    ss.kill(a**2)
    ss.kill(u**2)

    p_1 = ss.add_page({a: 0, u: 0})
    p_2 = ss.add_page({u: a})
    p_3 = ss.add_page()

    return {
        "name": "Hopf fibration S^1 -> S^3 -> S^2",
        "field": "GF(2)",
        "generators": {"a": a, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3},
        "ss": ss,
        "first_nonzero_differential": "d_2(u) = a",
        "expected_stable_page": "p_3",
    }


def build_example_3_universal_circle_bundle() -> Dict[str, object]:
    """Build the spectral sequence for S^1 -> ES^1 -> CP^infty."""
    c, u = symbols("c u")
    ss = SpectralSequence(
        GF(2),
        [c, u],
        [[2, 0], [0, 1]],
        [[1, 0], [-1, 1]],
    )
    ss.kill(u**2)

    p_1 = ss.add_page({c: 0, u: 0})
    p_2 = ss.add_page({u: c})
    p_3 = ss.add_page()

    return {
        "name": "Universal circle bundle S^1 -> ES^1 -> CP^infty",
        "field": "GF(2)",
        "generators": {"c": c, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3},
        "ss": ss,
        "first_nonzero_differential": "d_2(u) = c",
        "expected_stable_page": "p_3",
    }


def build_example_4_quaternionic_hopf() -> Dict[str, object]:
    """Build the spectral sequence for the quaternionic Hopf fibration."""
    a, u = symbols("a u")
    ss = SpectralSequence(
        GF(2),
        [a, u],
        [[4, 0], [0, 3]],
        [[1, 0], [-1, 1]],
    )
    ss.kill(a**2)
    ss.kill(u**2)

    p_1 = ss.add_page({a: 0, u: 0})
    p_2 = ss.add_page({a: 0, u: 0})
    p_3 = ss.add_page({a: 0, u: 0})
    p_4 = ss.add_page({u: a})
    p_5 = ss.add_page()

    return {
        "name": "Quaternionic Hopf fibration S^3 -> S^7 -> S^4",
        "field": "GF(2)",
        "generators": {"a": a, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3, "p_4": p_4, "p_5": p_5},
        "ss": ss,
        "first_nonzero_differential": "d_4(u) = a",
        "expected_stable_page": "p_5",
    }


def build_example_5_path_space_kz3() -> Dict[str, object]:
    """Build the mod-3 path-space example K(Z,2) -> PK(Z,3) -> K(Z,3)."""
    a, t = symbols("a t")
    ss = SpectralSequence(
        GF(3),
        [a, t],
        [[3, 0], [0, 2]],
        [[1, 0], [-1, 1]],
    )
    ss.kill(a**2)

    p_1 = ss.add_page({a: 0, t: 0})
    p_2 = ss.add_page({a: 0, t: 0})
    p_3 = ss.add_page({t: a})
    p_4 = ss.add_page()

    return {
        "name": "Path-space example K(Z,2) -> PK(Z,3) -> K(Z,3)",
        "field": "GF(3)",
        "generators": {"a": a, "t": t},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3, "p_4": p_4},
        "ss": ss,
        "first_nonzero_differential": "d_3(t) = a",
        "expected_stable_page": "p_4",
    }


def all_examples() -> List[ExampleRecord]:
    """Return all example records in a fixed order."""
    return [
        ExampleRecord(
            name="example_1_trivial_product",
            description="Trivial product fibration S^1 -> S^2 x S^1 -> S^2 over F_2.",
            build=build_example_1_trivial_product,
        ),
        ExampleRecord(
            name="example_2_hopf_fibration",
            description="Hopf fibration S^1 -> S^3 -> S^2 over F_2.",
            build=build_example_2_hopf_fibration,
        ),
        ExampleRecord(
            name="example_3_universal_circle_bundle",
            description="Universal circle bundle S^1 -> ES^1 -> CP^infty over F_2.",
            build=build_example_3_universal_circle_bundle,
        ),
        ExampleRecord(
            name="example_4_quaternionic_hopf",
            description="Quaternionic Hopf fibration S^3 -> S^7 -> S^4 over F_2.",
            build=build_example_4_quaternionic_hopf,
        ),
        ExampleRecord(
            name="example_5_path_space_kz3",
            description="Path-space model K(Z,2) -> PK(Z,3) -> K(Z,3) over F_3.",
            build=build_example_5_path_space_kz3,
        ),
    ]


def build_all_examples() -> Dict[str, Dict[str, object]]:
    """Construct and return all examples keyed by registry name."""
    output: Dict[str, Dict[str, object]] = {}
    for record in all_examples():
        output[record.name] = record.build()
    return output


if __name__ == "__main__":
    examples = build_all_examples()
    for key, data in examples.items():
        print(f"{key}: {data['name']}")
        print(f"  field: {data['field']}")
        print(f"  first nonzero differential: {data['first_nonzero_differential']}")
        print(f"  expected stable page: {data['expected_stable_page']}")
        print()
