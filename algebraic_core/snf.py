"""snf.py

Minimal API surface used by `spectral_sequence.py`.

This file intentionally contains only stubs (no implementation).
"""

from __future__ import annotations
from matrices import *
from sympy.matrices.normalforms import smith_normal_decomp

__all__ = ["SNF", "SNFMatrix"]
_verify = True


class SNF:
    @staticmethod
    def decomp(M: DMatrix):
        res = smith_normal_decomp(M.to_Matrix(), domain=M.domain)
        return [DMatrix.from_Matrix(_, domain=M.domain) for _ in res]

    @staticmethod
    def diag_divide(A: DMatrix, D: DMatrix) -> (DMatrix, list[int]):
        """
        Divide each row of A by corresponding diagonal element of D.
        Raise sympy.polys.polyerrors.ExactQuotientFailed when division fails.

        Return: One solution and the free rows
        """
        from sympy.polys.polyerrors import ExactQuotientFailed
        A, D = A.unify(D)
        domain = A.domain
        assert A.shape[0] == D.shape[0]
        rows = A.to_list()
        new_rows = []
        free_rows = []
        zero_row = [domain.zero] * A.shape[1]
        flag = False
        for i, row in enumerate(rows):
            if i < D.shape[1]:
                d = D[i, i].element
            else:
                # In the case that columns < rows, we need to stop early to ensure dimension compatibility
                # and also avoid index out of range.
                flag = True
                break

            if d != domain.zero:
                new_row = [domain.exquo(a, d) for a in row]
            else:
                new_row = zero_row
                for a in row:
                    if a != domain.zero:
                        raise ExactQuotientFailed(a, d, dom=domain)
                free_rows.append(i)
            new_rows.append(new_row)

        # In the case that rows < columns, we need to pad additional zero rows to ensure dimension compatibility.
        new_rows.extend([zero_row] * (D.shape[1] - len(new_rows)))

        res = DMatrix.from_list(new_rows, domain=domain)
        if D * res != A:
            if flag:
                raise ExactQuotientFailed(A, D)
            else:
                raise AssertionError(A, D, res, D*res)
        return res, free_rows

    @staticmethod
    def solve(T: DMatrix, A: DMatrix, U=None, D=None, V=None) -> tuple[DMatrix, DMatrix] | None:
        """
        Solve T = AX for X.

        Return: one solution and columns of V whose span gives the kernel.
        Note: the free_rows returned by diag_divide are where D=A.D has zero diagonal elements,
        so it only depends on A and is exactly the kernel of A.

        """

        from sympy.polys.polyerrors import ExactQuotientFailed
        try:
            if U is None or D is None or V is None:
                D, U, V = SNF.decomp(A)
            S, free_rows = SNF.diag_divide(U * T, D)
            res = V * S
            assert A * res == T
            return res, V.extract_columns(free_rows)
        except ExactQuotientFailed:
            return None

    @staticmethod
    def align(A: DMatrix, B: DMatrix, _X=None):
        """
        Return P, Q, D such that AP = BQD, where P, Q are invertible and D is diagonal.
        It requires that A=BX is solvable.
        """
        if _X is None:
            _X = SNF.solve(A, B)[0]
        if _X is not None:
            D, U, V = SNF.decomp(_X)
            assert A * V == B * U.inv_den()[0] * D
            return V, U.inv_den()[0], D
        assert False

    @staticmethod
    def kernel_of(A: DMatrix):
        D, U, V = SNF.decomp(A)
        free_columns = [i for i in range(D.shape[0]) if D[i, i].element == A.domain.zero]
        return V.extract_columns(free_columns)

    @staticmethod
    def spans(M: DMatrix, V: DMatrix):
        assert V.shape == (M.shape[0], 1)
        return SNF.solve(V, M) is not None


class SNFMatrix(DMatrix):
    """Matrix that caches SNF (or field rank/nullspace) for kernel/span queries."""
    @classmethod
    def from_rep(cls, rep) -> SNFMatrix:
        """The __new__ in domain matrix invokes from_rep, so we only need to update this."""
        instance = super().from_rep(rep)
        assert isinstance(instance, SNFMatrix)

        domain = rep.domain
        content = IMatrix(rep.to_list())
        if _verify:
            assert domain.is_PID
        D, U, V = smith_normal_decomp(content, domain=domain)
        instance.D = DMatrix.from_Matrix(D, domain=domain)
        instance.U = DMatrix.from_Matrix(U, domain=domain)
        instance.V = DMatrix.from_Matrix(V, domain=domain)

        return instance

    def solve(self, T: DMatrix) -> DMatrix | None:
        """
        Solve T = self * X for X. Return None if not solvable
        """
        return SNF.solve(T, self, self.U, self.D, self.V)

    def kernel(self):
        return SNF.kernel_of(self)

    def spans(self, v: DMatrix):
        assert v.shape == (self.shape[0], 1)
        return self.solve(v) is not None

    def align(self, M):
        return SNF.align(M, self, _X=self.solve(M))

    def __str__(self):
        from sympy import MatrixBase
        return "SNF" + MatrixBase.__str__(self.to_Matrix())
