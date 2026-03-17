from __future__ import annotations

from src.matrices import *

__all__ = ["SNF", "SNFMatrix"]
_verify = True


class SNF:
    @staticmethod
    def _identity_dense(size: int, domain):
        one = domain.one
        zero = domain.zero
        return [[one if i == j else zero for j in range(size)] for i in range(size)]

    @staticmethod
    def _swap_rows_dense(M, i: int, j: int):
        M[i], M[j] = M[j], M[i]

    @staticmethod
    def _swap_cols_dense(M, i: int, j: int):
        if i == j:
            return
        for r in range(len(M)):
            M[r][i], M[r][j] = M[r][j], M[r][i]

    @staticmethod
    def _combine_rows_dense(M, i: int, j: int, a, b, c, d):
        row_i = M[i][:]
        row_j = M[j][:]
        M[i] = [a * x + b * y for x, y in zip(row_i, row_j)]
        M[j] = [c * x + d * y for x, y in zip(row_i, row_j)]

    @staticmethod
    def _combine_cols_dense(M, i: int, j: int, a, b, c, d):
        col_i = [row[i] for row in M]
        col_j = [row[j] for row in M]
        for r in range(len(M)):
            M[r][i] = a * col_i[r] + b * col_j[r]
            M[r][j] = c * col_i[r] + d * col_j[r]

    @staticmethod
    def _find_nonzero_in_block(A, k: int, domain):
        zero = domain.zero
        m = len(A)
        n = len(A[0]) if m else 0
        for i in range(k, m):
            for j in range(k, n):
                if A[i][j] != zero:
                    return i, j
        return None

    @staticmethod
    def _unit_inverse(domain, x):
        if not domain.is_unit(x):
            raise ValueError(f"Expected a unit in domain {domain}, got {x}.")
        return domain.exquo(domain.one, x)

    @staticmethod
    def invert_unimodular(M: DMatrix) -> DMatrix:
        """
        Return the exact inverse of a unimodular matrix over the same domain.

        DomainMatrix.inv_den() gives (N, d) with N * M = d * I.
        For unimodular M, d is a unit, so M^{-1} = d^{-1} * N in-domain.
        """
        inv_num, den = M.inv_den()
        den_inv = SNF._unit_inverse(M.domain, den)
        return inv_num * den_inv

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
        """
        Smith normal decomposition over a PID using unimodular row/column operations.

        Return D, U, V such that U * M * V = D.
        """
        m, n = M.shape
        domain = M.domain
        zero = domain.zero
        one = domain.one

        if m == 0 or n == 0:
            D = DMatrix.zeros((m, n), domain)
            U = DMatrix.eye((m, m), domain)
            V = DMatrix.eye((n, n), domain)
            return [D, U, V]

        if all(x == zero for x in M.to_list_flat()):
            D = DMatrix.zeros((m, n), domain)
            U = DMatrix.eye((m, m), domain)
            V = DMatrix.eye((n, n), domain)
            return [D, U, V]

        A = [row[:] for row in M.to_list()]
        U = SNF._identity_dense(m, domain)
        V = SNF._identity_dense(n, domain)

        k = 0
        while k < min(m, n):
            pivot_pos = SNF._find_nonzero_in_block(A, k, domain)
            if pivot_pos is None:
                break
            i0, j0 = pivot_pos
            if i0 != k:
                SNF._swap_rows_dense(A, k, i0)
                SNF._swap_rows_dense(U, k, i0)
            if j0 != k:
                SNF._swap_cols_dense(A, k, j0)
                SNF._swap_cols_dense(V, k, j0)

            guard = 0
            while True:
                guard += 1
                if guard > (m + n + 5) * (m + n + 5):
                    raise RuntimeError("SNF decomposition did not converge; possible non-Euclidean behavior.")

                changed = False
                pivot = A[k][k]

                if pivot == zero:
                    nxt = SNF._find_nonzero_in_block(A, k, domain)
                    if nxt is None:
                        break
                    i1, j1 = nxt
                    if i1 != k:
                        SNF._swap_rows_dense(A, k, i1)
                        SNF._swap_rows_dense(U, k, i1)
                    if j1 != k:
                        SNF._swap_cols_dense(A, k, j1)
                        SNF._swap_cols_dense(V, k, j1)
                    changed = True
                    if changed:
                        continue

                # Clear column k below pivot.
                for i in range(k + 1, m):
                    b = A[i][k]
                    if b == zero:
                        continue
                    s, t, g = domain.gcdex(pivot, b)
                    u = domain.exquo(pivot, g)
                    v = domain.exquo(b, g)
                    # [[s, t], [-v, u]] is unimodular because s*u + t*v = 1
                    SNF._combine_rows_dense(A, k, i, s, t, -v, u)
                    SNF._combine_rows_dense(U, k, i, s, t, -v, u)
                    pivot = A[k][k]
                    changed = True

                # Clear row k to the right of pivot.
                for j in range(k + 1, n):
                    b = A[k][j]
                    if b == zero:
                        continue
                    s, t, g = domain.gcdex(pivot, b)
                    u = domain.exquo(pivot, g)
                    v = domain.exquo(b, g)
                    SNF._combine_cols_dense(A, k, j, s, t, -v, u)
                    SNF._combine_cols_dense(V, k, j, s, t, -v, u)
                    pivot = A[k][k]
                    changed = True

                pivot = A[k][k]
                if pivot == zero:
                    continue

                # Ensure pivot divides all entries in the lower-right block.
                violating = None
                if not domain.is_unit(pivot):
                    for i in range(k + 1, m):
                        for j in range(k + 1, n):
                            if domain.rem(A[i][j], pivot) != zero:
                                violating = (i, j)
                                break
                        if violating is not None:
                            break

                if violating is not None:
                    i, _ = violating
                    # Introduce a violating entry into pivot row, then re-reduce.
                    # This unimodular operation keeps rank/invariants and decreases pivot up to associates.
                    SNF._combine_rows_dense(A, k, i, one, one, zero, one)
                    SNF._combine_rows_dense(U, k, i, one, one, zero, one)
                    changed = True

                if not changed:
                    break

            if domain.is_ZZ and A[k][k] < 0:
                A[k] = [-x for x in A[k]]
                U[k] = [-x for x in U[k]]

            k += 1

        D = DMatrix.from_list(A, domain=domain)
        U = DMatrix.from_list(U, domain=domain)
        V = DMatrix.from_list(V, domain=domain)

        assert U * M * V == D
        return [D, U, V]

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
                raise AssertionError(A, D, res, D * res)
        return res, free_columns

    @staticmethod
    def solve(T: DMatrix, A: DMatrix, U=None, D=None, V=None) -> tuple[DMatrix, DMatrix] | None:
        """
        Solve T = AX for X.

        Return: one solution and columns of V whose span gives the kernel.
        Note: free_columns from diag_divide are where D has zero diagonal elements,
        so it depends only on A and describes the kernel directions.
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
            solved = SNF.solve(A, B)
            if solved is None:
                raise ValueError("Cannot align: A is not in the column span of B.")
            _X = solved[0]
        if A != B * _X:
            raise ValueError("Cannot align: expected A = B * X.")

        D, U, V = SNF.decomp(_X)
        U_inv = SNF.invert_unimodular(U)
        assert _X * V == U_inv * D
        assert A * V == B * U_inv * D
        return V, U_inv, D

    @staticmethod
    def kernel_of(A: DMatrix):
        D, U, V = SNF.decomp(A)
        free_columns = SNF._kernel_free_columns(D)
        return V.extract_columns(free_columns)

    @staticmethod
    def spans(M: DMatrix, V: DMatrix):
        assert V.shape == (M.shape[0], 1)
        return SNF.solve(V, M) is not None


class SNFMatrix(DMatrix):
    """Matrix that caches SNF decomposition for kernel/span queries."""

    @classmethod
    def from_rep(cls, rep) -> SNFMatrix:
        """The __new__ in domain matrix invokes from_rep, so we only need to update this."""
        instance = super().from_rep(rep)
        assert isinstance(instance, SNFMatrix)

        domain = rep.domain
        if _verify:
            assert domain.is_PID

        M = DMatrix.from_Matrix(IMatrix(rep.to_list()), domain=domain)
        D, U, V = SNF.decomp(M)
        instance.D = D
        instance.U = U
        instance.V = V
        return instance

    def solve(self, T: DMatrix) -> DMatrix | None:
        """
        Solve T = self * X for X. Return None if not solvable.
        """
        return SNF.solve(T, self, self.U, self.D, self.V)

    def kernel(self):
        return SNF.kernel_of(self)

    def spans(self, v: DMatrix):
        assert v.shape == (self.shape[0], 1)
        return self.solve(v) is not None

    def align(self, M):
        solved = self.solve(M)
        if solved is None:
            raise ValueError("Cannot align: matrix is not in span.")
        return SNF.align(M, self, _X=solved[0])

    def __str__(self):
        from sympy import MatrixBase

        return "SNF" + MatrixBase.__str__(self.to_Matrix())
