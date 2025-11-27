import unittest
from sympy.polys.domains import ZZ, QQ
from sympy import Rational

from AlgebraicCore import SNFMatrix, divides, Module, Differential
from utilities import Matrix


class TestAlgebraicCore(unittest.TestCase):
    def test_divides_over_zz_and_field(self):
        # PID behavior
        self.assertTrue(divides(6, 3, ZZ))
        self.assertFalse(divides(5, 3, ZZ))
        self.assertTrue(divides(0, 0, ZZ))
        self.assertFalse(divides(5, 0, ZZ))
        # Field fast path: any nonzero divisor divides anything
        self.assertTrue(divides(Rational(1, 2), Rational(3, 7), QQ))
        self.assertTrue(divides(0, 2, QQ))  # 2 divides 0

    def test_snfmatrix_kernel_and_span_over_zz(self):
        A = SNFMatrix([[1, 1], [1, 1]], domain=ZZ)
        self.assertEqual(A.r, 1)
        K = A.get_kernel()
        self.assertEqual(K.shape, (2, 1))
        # A * kernel == 0
        self.assertEqual(A * K, Matrix.zeros(A.rows, K.cols))
        # Column span checks
        self.assertTrue(A.column_spans(Matrix([[2], [2]])))
        self.assertFalse(A.column_spans(Matrix([[1], [0]])))

    def test_snfmatrix_field_path_rank_and_span(self):
        A = SNFMatrix([[1, 1], [1, 1]], domain=QQ)
        self.assertEqual(A.r, 1)
        # Kernel dimension 1 over QQ
        K = A.get_kernel()
        self.assertEqual(K.shape[1], 1)
        self.assertEqual(A * K, Matrix.zeros(A.rows, K.cols))
        # Span checks via rank test
        self.assertTrue(A.column_spans(Matrix([[2], [2]])))
        self.assertFalse(A.column_spans(Matrix([[1], [0]])))

    def test_module_verify_and_classify(self):
        span_set = [Matrix([[1], [2]]), Matrix([[1], [2]])]
        zero_set = [Matrix([[2], [4]])]
        m = Module(ZZ, 2, span_set, zero_set)

        self.assertEqual(m.classify(Matrix([[0], [0]])), "in_zero_span")
        self.assertEqual(m.classify(Matrix([[2], [4]])), "in_zero_span")
        self.assertEqual(m.classify(Matrix([[1], [2]])), "in_span_not_zero")
        self.assertEqual(m.classify(Matrix([[1], [0]])), "outside_span")

    def test_module_empty_zero_set_defaults_to_zero(self):
        span_set = [Matrix([[1], [0]]), Matrix([[0], [1]])]
        zero_set: list[Matrix] = []
        m = Module(ZZ, 2, span_set, zero_set)

        self.assertEqual(m.classify(Matrix([[0], [0]])), "in_zero_span")
        self.assertEqual(m.classify(Matrix([[1], [0]])), "in_span_not_zero")
        self.assertEqual(m.classify(Matrix([[1], [1]])), "in_span_not_zero")
        self.assertEqual(m.classify(Matrix([[2], [3]])), "in_span_not_zero")
        self.assertEqual(m.classify(Matrix([[5], [0]])), "in_span_not_zero")

    def test_differential_builds_matrix_and_prunes(self):
        # from_module: span {e1, e2}, with K = span{e2}
        e1 = Matrix([[1], [0]])
        e2 = Matrix([[0], [1]])
        from_mod = Module(ZZ, 2, [e1, e2], [e2])
        to_mod = Module(ZZ, 1, [Matrix([[1]])], [])
        non_zero = {e1: Matrix([[3]])}  # d(e1) = 3, d(e2)=0
        d = Differential(from_mod, to_mod, non_zero)
        self.assertTrue(d.compute_matrix())
        self.assertEqual(d.matrix.shape, (1, 2))
        self.assertEqual(d.matrix, Matrix([[3, 0]]))
        # non_zero_info should be pruned to basis generators (e1 only)
        self.assertEqual(len(d.non_zero_info), 1)
        self.assertIn(e1, d.non_zero_info)

    def test_differential_insufficient_generators_returns_false(self):
        e1 = Matrix([[1], [0]])
        e2 = Matrix([[0], [1]])
        from_mod = Module(ZZ, 2, [e1, e2], [])
        to_mod = Module(ZZ, 1, [Matrix([[1]])], [])
        non_zero = {e1: Matrix([[5]])}  # Only e1 given; e2 missing
        d = Differential(from_mod, to_mod, non_zero)
        self.assertFalse(d.compute_matrix())
        self.assertIsNone(d.matrix)


if __name__ == "__main__":
    unittest.main()
