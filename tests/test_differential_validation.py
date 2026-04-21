from __future__ import annotations

import sys
from pathlib import Path

import pytest
from sympy import GF
from sympy.abc import a, u


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for path in (SRC, ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from src.spectral_sequence import SpectralSequence  # noqa: E402


def test_rejects_nonzero_image_of_relation_trivial_source():
    ss = SpectralSequence(GF(2), [a, u], [[2, 0], [0, 1]], [[1, 0], [-1, 1]])
    ss.kill(a**2)
    ss.kill(u**2)
    ss.add_page({a: 0, u: 0})

    with pytest.raises(ValueError, match="must be zero"):
        ss.add_page({u**2: a * u})


def test_accepts_zero_image_of_relation_trivial_source():
    ss = SpectralSequence(GF(2), [a, u], [[2, 0], [0, 1]], [[1, 0], [-1, 1]])
    ss.kill(a**2)
    ss.kill(u**2)
    ss.add_page({a: 0, u: 0})

    page = ss.add_page({u**2: 0})
    assert page.page_num == 2
