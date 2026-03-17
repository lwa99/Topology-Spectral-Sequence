from __future__ import annotations
from src.matrices import *
from sympy.matrices.normalforms import smith_normal_decomp

__all__ = ["SNF", "SNFMatrix"]
_verify = True


class SNF:
    @staticmethod
    def _kernel_free_columns(D: DMatrix):
        """
        Return indices of free columns in the diagonal Smith form matrix D.

        For D of shape (m, n), variable i is free iff:
        - i >= m (extra columns beyond number of equations), or
        - i < min(m, n) and D[i, i] == 0.
        """
        m, n = D.shape
        diag_len = min(m, n)
        free_columns = []
        for i in range(n):
            if i >= m:
                free_columns.append(i)
            elif i < diag_len and D[i, i].element == D.domain.zero:
                free_columns.append(i)
        return free_columns

    @staticmethod
    def decomp(M: DMatrix):
        from sympy.polys.polyerrors import NotInvertible
        m, n = M.shape
        domain = M.domain

        if m == 0 or n == 0:
            D = DMatrix.zeros((m, n), domain)
            U = DMatrix.eye((m, m), domain)
            V = DMatrix.eye((n, n), domain)
            return [D, U, V]

        if all(x == domain.zero for x in M.to_list_flat()):
            D = DMatrix.zeros((m, n), domain)
            U = DMatrix.eye((m, m), domain)
            V = DMatrix.eye((n, n), domain)
            return [D, U, V]

        try:
            res = smith_normal_decomp(M.to_Matrix(), domain=domain)
            return [DMatrix.from_Matrix(_, domain=domain) for _ in res]
        except NotInvertible:
            # SymPy may fail on some degenerate inputs over fields; the all-zero case was handled above.
            raise

    @staticmethod
    def diag_divide(A: DMatrix, D: DMatrix) -> (DMatrix, list[int]):
        """
        Divide each row of A by corresponding diagonal element of D.
        Raise sympy.polys.polyerrors.ExactQuotientFailed when division fails.

        Return: one solution and free variable indices (columns in the unknown matrix).
        """
        from sympy.polys.polyerrors import ExactQuotientFailed
        A, D = A.unify(D)
        domain = A.domain
        assert A.shape[0] == D.shape[0]
        if D.shape[1] == 0:
            for row in A.to_list():
                for a in row:
                    if a != domain.zero:
                        raise ExactQuotientFailed(A, D)
            return DMatrix.zeros((0, A.shape[1]), domain), []

        rows = A.to_list()
        new_rows = []
        free_columns = []
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
                free_columns.append(i)
            new_rows.append(new_row)

        # In the case that rows < columns, we need to pad additional zero rows to ensure dimension compatibility.
        new_rows.extend([zero_row] * (D.shape[1] - len(new_rows)))
        if D.shape[1] > D.shape[0]:
            free_columns.extend(range(D.shape[0], D.shape[1]))

        res = DMatrix.from_list(new_rows, domain=domain)
        if D * res != A:
            if flag:
                raise ExactQuotientFailed(A, D)
            else:
                raise AssertionError(A, D, res, D*res)
        return res, free_columns

    @staticmethod
    def solve(T: DMatrix, A: DMatrix, U=None, D=None, V=None) -> tuple[DMatrix, DMatrix] | None:
        """
        Solve T = AX for X.

        Return: one solution and columns of V whose span gives the kernel.
        Note: the free_rows returned by diag_divide are where D=A.D has zero diagonal elements,
        so it only depends on A and is exactly the kernel of A.

        """

        from sympy.polys.polyerrors import ExactQuotientFailed, NotInvertible
        T, A = T.unify(A)
        domain = A.domain

        # Zero-row system: every RHS is zero in that ambient shape, so solvable.
        # One canonical solution is X = 0; kernel is all columns in the source.
        if A.shape[0] == 0:
            assert T.shape[0] == 0
            sol = DMatrix.zeros((A.shape[1], T.shape[1]), domain)
            ker = DMatrix.eye((A.shape[1], A.shape[1]), domain)
            return sol, ker

        # Zero-column system: AX is always zero; solvable iff RHS is zero.
        if A.shape[1] == 0:
            if any(x != domain.zero for x in T.to_list_flat()):
                return None
            sol = DMatrix.zeros((0, T.shape[1]), domain)
            ker = DMatrix.zeros((0, 0), domain)
            return sol, ker

        try:
            if U is None or D is None or V is None:
                D, U, V = SNF.decomp(A)
            S, free_columns = SNF.diag_divide(U * T, D)
            assert A * V * S == T, (A * V * S, T)
            res = V * S
            return res, V.extract_columns(free_columns)
        except (ExactQuotientFailed, NotInvertible):
            return None

    @staticmethod
    def align(A: DMatrix, B: DMatrix, _X=None) -> tuple[DMatrix, DMatrix, DMatrix]:
        """
        Return P, Q, D such that AP = BQD, where P, Q are invertible and D is diagonal.
        It requires that A=BX is solvable.
        """
        if _X is None:
            _X = SNF.solve(A, B)[0]
        assert A == B * _X
        D, U, V = SNF.decomp(_X)
        print(U, U.inv_den())
        assert _X * V == U.inv_den()[0] * D
        assert A * V == B * U.inv_den()[0] * D
        return V, U.inv_den()[0], D

    @staticmethod
    def kernel_of(A: DMatrix):
        if A.domain.is_Field:
            # Over fields, nullspace via linear algebra is robust and avoids SNF edge cases in SymPy.
            basis = A.to_Matrix().nullspace()
            if len(basis) == 0:
                return DMatrix.zeros((A.shape[1], 0), A.domain)
            cols = [DMatrix.from_Matrix(v, domain=A.domain) for v in basis]
            return DMatrix.static_hstack(*cols)
        D, U, V = SNF.decomp(A)
        free_columns = SNF._kernel_free_columns(D)
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
        if _verify:
            assert domain.is_PID
        if domain.is_Field:
            # Over fields, avoid eager Smith decomposition: SymPy has edge-case failures on some inputs.
            instance.D = None
            instance.U = None
            instance.V = None
            return instance
        M = DMatrix.from_Matrix(IMatrix(rep.to_list()), domain=domain)
        D, U, V = SNF.decomp(M)
        instance.D = D
        instance.U = U
        instance.V = V

        return instance

    def solve(self, T: DMatrix) -> DMatrix | None:
        """
        Solve T = self * X for X. Return None if not solvable
        """
        if self.domain.is_Field:
            return SNF.solve(T, self)
        return SNF.solve(T, self, self.U, self.D, self.V)

    def kernel(self):
        return SNF.kernel_of(self)

    def spans(self, v: DMatrix):
        assert v.shape == (self.shape[0], 1)
        return self.solve(v) is not None

    def align(self, M):
        return SNF.align(M, self, _X=self.solve(M)[0])

    def __str__(self):
        from sympy import MatrixBase
        return "SNF" + MatrixBase.__str__(self.to_Matrix())
