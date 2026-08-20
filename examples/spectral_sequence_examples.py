from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from sympy import GF
from sympy.abc import symbols

from src.spectral_sequence import SpectralSequence


@dataclass(frozen=True)
class ExampleRecord:
    name: str
    build: Callable[[], dict]


def _new_ss(domain, generators, bidegrees):
    return SpectralSequence(
        domain,
        generators,
        bidegrees,
        [[1, 0], [-1, 1]],
    )


def example_1_trivial_product():
    a, b = symbols("a b")
    ss = _new_ss(GF(2), [a, b], [[2, 0], [0, 1]])
    ss.kill(a**2, b**2)
    p_1 = ss.add_page({a: 0, b: 0})
    p_2 = ss.add_page()
    return {
        "name": "Trivial product",
        "field": "GF(2)",
        "ss": ss,
        "generators": {"a": a, "b": b},
        "pages": {"p_1": p_1, "p_2": p_2},
        "expected_stable_page": "p_2",
    }


def example_2_hopf_fibration():
    a, u = symbols("a u")
    ss = _new_ss(GF(2), [a, u], [[2, 0], [0, 1]])
    ss.kill(a**2, u**2)
    p_1 = ss.add_page({a: 0, u: 0})
    p_2 = ss.add_page({u: a, a: 0})
    p_3 = ss.add_page()
    return {
        "name": "Hopf fibration",
        "field": "GF(2)",
        "ss": ss,
        "generators": {"a": a, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3},
        "expected_stable_page": "p_3",
    }


def example_3_universal_circle_bundle():
    c, u = symbols("c u")
    ss = _new_ss(GF(2), [c, u], [[2, 0], [0, 1]])
    ss.kill(u**2)
    p_1 = ss.add_page({c: 0, u: 0})
    p_2 = ss.add_page({u: c, c: 0})
    p_3 = ss.add_page()
    return {
        "name": "Universal circle bundle",
        "field": "GF(2)",
        "ss": ss,
        "generators": {"c": c, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3},
        "expected_stable_page": "p_3",
    }


def example_4_quaternionic_hopf():
    a, u = symbols("a u")
    ss = _new_ss(GF(2), [a, u], [[4, 0], [0, 3]])
    ss.kill(a**2, u**2)
    p_1 = ss.add_page({a: 0, u: 0})
    p_2 = ss.add_page({a: 0, u: 0})
    p_3 = ss.add_page({a: 0, u: 0})
    p_4 = ss.add_page({u: a, a: 0})
    p_5 = ss.add_page()
    return {
        "name": "Quaternionic Hopf fibration",
        "field": "GF(2)",
        "ss": ss,
        "generators": {"a": a, "u": u},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3, "p_4": p_4, "p_5": p_5},
        "expected_stable_page": "p_5",
    }


def example_5_path_space_kz3():
    a, t = symbols("a t")
    ss = _new_ss(GF(3), [a, t], [[3, 0], [0, 2]])
    ss.kill(a**2)
    p_1 = ss.add_page({a: 0, t: 0})
    p_2 = ss.add_page({a: 0, t: 0})
    p_3 = ss.add_page({t: a, a: 0})
    p_4 = ss.add_page()
    return {
        "name": "Path-space style example over GF(3)",
        "field": "GF(3)",
        "ss": ss,
        "generators": {"a": a, "t": t},
        "pages": {"p_1": p_1, "p_2": p_2, "p_3": p_3, "p_4": p_4},
        "expected_stable_page": "p_4",
    }


def all_examples():
    return [
        ExampleRecord("example_1_trivial_product", example_1_trivial_product),
        ExampleRecord("example_2_hopf_fibration", example_2_hopf_fibration),
        ExampleRecord("example_3_universal_circle_bundle", example_3_universal_circle_bundle),
        ExampleRecord("example_4_quaternionic_hopf", example_4_quaternionic_hopf),
        ExampleRecord("example_5_path_space_kz3", example_5_path_space_kz3),
    ]


def build_all_examples():
    return {record.name: record.build() for record in all_examples()}


if __name__ == "__main__":
    built = build_all_examples()
    for key, data in built.items():
        print(f"{key}: {data['name']} ({data['field']}) -> {data['expected_stable_page']}")
