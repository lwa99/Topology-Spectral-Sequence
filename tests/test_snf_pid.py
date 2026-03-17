from __future__ import annotations

import sys
from pathlib import Path

import pytest
from sympy import GF, ZZ


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for path in (SRC, ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from matrices import DMatrix  # noqa: E402
from snf import SNF  # noqa: E402


@pytest.mark.parametrize(
    "domain,rows",
    [
        (ZZ, [[2, 4], [6, 8]]),
        (GF(3), [[0, 0, 0], [0, 0, 0], [0, 0, -1]]),
    ],
)
def test_align_handles_unit_denominator(domain, rows):
    A = DMatrix.from_list(rows, domain)
    n = A.shape[0]
    B = DMatrix.eye((n, n), domain)

    P, Q, D = SNF.align(A, B)
    Q_inv = SNF.invert_unimodular(Q)

    assert A * P == B * Q * D
    assert Q * Q_inv == DMatrix.eye((Q.shape[0], Q.shape[0]), domain)
    assert Q_inv * Q == DMatrix.eye((Q.shape[1], Q.shape[1]), domain)


@pytest.mark.parametrize(
    "rows",
    [
        [[0, 0, 0, 1]],
        [[0, 0, 0, 2]],
    ],
)
def test_custom_decomp_handles_gf_corner_case(rows):
    """
    SymPy smith_normal_decomp is unstable on these matrices over GF(3).
    Our custom decomposition should still produce U*M*V = D.
    """
    domain = GF(3)
    M = DMatrix.from_list(rows, domain)
    D, U, V = SNF.decomp(M)

    assert U * M * V == D
    assert D.shape == M.shape
    # Rank-1 row vector over a field should have one nonzero diagonal entry.
    diag = D.diagonal()
    assert len([x for x in diag if x != domain.zero]) == 1
