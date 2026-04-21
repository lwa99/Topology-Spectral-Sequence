from __future__ import annotations

import sys
from pathlib import Path

from sympy import ZZ
from sympy.abc import a, t


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for path in (SRC, ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from src.matrices import DV  # noqa: E402
from src.page_and_module import Module  # noqa: E402
from src.spectral_sequence import SpectralSequence  # noqa: E402


def test_structural_information_uses_basis_not_raw_spanning_set():
    ss = SpectralSequence(ZZ, [a], [[1], [0]], [[1, 0], [-1, 1]])
    p1 = ss.add_page({a: 0})

    module = Module(
        p1,
        p1._normalize_bidegree((1, 0)),
        [DV([1], ZZ), DV([1], ZZ)],
        [],
    )

    info = module.get_structural_information()
    assert info is not None
    gens, torsion = info

    assert len(gens) == 1
    assert torsion == [ZZ.zero]
    assert module.classify(gens[0]) == 1


def test_structural_information_omits_unit_torsion_zero_summands():
    ss = SpectralSequence(ZZ, [a, t], [[3, 0], [0, 2]], [[1, 0], [-1, 1]])
    ss.kill(a**2)
    ss.add_page({a: 0, t: 0})
    p2 = ss.add_page({a: 0, t: 0})

    info = p2[6, 0].get_structural_information()
    assert info == ([], [])


def test_structural_information_keeps_actual_torsion_summands():
    ss = SpectralSequence(ZZ, [a, t], [[3, 0], [0, 2]], [[1, 0], [-1, 1]])
    ss.kill(a**2)
    ss.add_page({a: 0, t: 0})
    ss.add_page({a: 0, t: 0})
    ss.add_page({t: a, 2 * a: 0})
    p4 = ss.add_page()

    info = p4[3, 2].get_structural_information()
    assert info is not None
    gens, torsion = info

    assert len(gens) == 1
    assert torsion == [ZZ(2)]
    assert p4[3, 2].classify(gens[0]) == 1
