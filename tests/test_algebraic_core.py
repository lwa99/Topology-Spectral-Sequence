import pytest
from sympy.polys.domains import ZZ

from AlgebraicCore import SNFMatrix, divides, Module, _hstack
from utilities import Matrix
from sympy import Rational
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.domains import GF
from sympy import Matrix as SMatrix


def test_divides_basic_zz():
    assert divides(6, 3, ZZ) is True
    assert divides(5, 3, ZZ) is False
    assert divides(-8, 4, ZZ) is True
    assert divides(8, -4, ZZ) is True


def test_divides_zero_divisor_behavior():
    # Only 0 divides 0
    assert divides(0, 0, ZZ) is True
    assert divides(5, 0, ZZ) is False


def test_snfmatrix_decomposition_and_rank():
    a = SNFMatrix([[2, 4], [6, 8]], domain=ZZ)
    # Validate U * A * V == D
    assert a._U * a * a._V == a._D
    # Rank equals number of nonzero diagonal entries of D
    nonzero_diagonals = sum(1 for i in range(min(a.rows, a.cols)) if a._D[i, i] != 0)
    assert a.r == nonzero_diagonals
    # Repr contains the data for quick sanity
    assert "SNFMatrix" in repr(a)


def test_column_spans_on_diagonal_matrix():
    # A = diag(2, 4) over ZZ; span is {(2a, 4b)}
    A = SNFMatrix([[2, 0], [0, 4]], domain=ZZ)
    assert A.column_spans(Matrix([[6], [8]])) is True
    assert A.column_spans(Matrix([[0], [0]])) is True
    assert A.column_spans(Matrix([[2], [2]])) is False


def test_column_spans_single_column():
    # A spans multiples of column [1, 0]^T
    A = SNFMatrix([[1], [0]], domain=ZZ)
    assert A.column_spans(Matrix([[5], [0]])) is True
    assert A.column_spans(Matrix([[1], [1]])) is False


def test_module_zero_set_validation_passes():
    # span_set columns generate multiples of [1, 2]^T, zero_set contains [1, 4]^T = 2 * [1, 2]^T
    span_set = [Matrix([[1], [2]]), Matrix([[1], [2]])]
    zero_set = [Matrix([[2], [4]])]
    Module(ZZ, 2, span_set, zero_set)  # should not raise


def test_module_zero_set_validation_fails():
    span_set = [Matrix([[1], [0]])]
    zero_set = [Matrix([[0], [1]])]  # not in the ZZ-span of [1, 0]^T
    with pytest.raises(AssertionError):
        Module(ZZ, 2, span_set, zero_set)


def test_column_spans_rank_deficient_lower_block_violation():
    # A has rank 1; the lower block condition must force b'_{2} = 0.
    A = SNFMatrix([[1, 0], [0, 0]], domain=ZZ)
    e2 = Matrix([[0], [1]])
    # Choose b so that U*b = e2 (violates lower block zero condition)
    invU = A._U.inv()
    b = invU * e2
    assert A.column_spans(b) is False


def test_hstack_helper_matches_manual():
    c1 = Matrix([[1], [2]])
    c2 = Matrix([[3], [4]])
    H1 = _hstack(ZZ, c1, c2)
    H2 = SNFMatrix(Matrix.hstack(c1, c2), domain=ZZ)
    assert H1 == H2


def test_uncoercible_entry_raises():
    # Using a Rational that is not in ZZ should cause a coercion failure
    with pytest.raises(CoercionFailed):
        SNFMatrix([[Rational(1, 2)]], domain=ZZ)


def test_divides_in_finite_field():
    F5 = GF(5)
    # In a field, every nonzero divides every nonzero; 0 divides only 0
    assert divides(3, 2, F5) is True
    assert divides(1, 4, F5) is True
    assert divides(0, 1, F5) is True  # 1 divides 0
    assert divides(2, 0, F5) is False


def test_column_spans_random_constructions_are_detected():
    # Randomized property: if b is constructed as A * x for integer x, column_spans must be True
    import random
    random.seed(0)
    for _ in range(10):
        rows = 2
        cols = 3
        data = [[random.randint(-2, 2) for _ in range(cols)] for _ in range(rows)]
        A = SNFMatrix(data, domain=ZZ)
        x = Matrix([[random.randint(-2, 2)] for _ in range(cols)])
        b = A * x
        assert A.column_spans(b) is True


def test_get_kernel_dimension_matches_nullspace():
    # A few matrices with varying kernel dimensions
    cases = [
        [[1, 2], [2, 4]],      # rank 1, kernel dim 1
        [[1, 2, 3]],           # rank 1, kernel dim 2
        [[2, 0], [0, 3]],      # full rank, kernel dim 0
        [[1, 1, 1], [1, 1, 1]] # rank 1, kernel dim 2
    ]
    for data in cases:
        A = SNFMatrix(data)
        K = A.get_kernel()
        # Nullspace via SymPy (over Q)
        NS_list = SMatrix(data).nullspace()
        null_cols = len(NS_list)
        assert K.cols == null_cols
        # Check A * K == 0
        if K.cols == 0:
            # A * K should be empty with correct shape
            assert K.rows == A.cols
        else:
            Z = A * K
            assert Z == Z*0  # quick structural equality check
            for i in range(Z.rows):
                for j in range(Z.cols):
                    assert Z[i, j] == 0


def test_get_kernel_entries_coerce_to_domain_ZZ():
    cases = [
        [[1, 2], [2, 4]],
        [[1, 2, 3]],
        [[0, 0], [0, 0]],
    ]
    for data in cases:
        A = SNFMatrix(data)
        K = A.get_kernel()
        for entry in K:
            # Should be coercible into ZZ
            A.domain.convert(entry)
