from sympy import ImmutableMatrix as IMatrix, Matrix, MatrixBase
from sympy.polys.matrices import DomainMatrix as _DMatrix
from collections.abc import Iterable

__all__ = ["IM", "IV", "DM", "DV", "IMatrix", "DMatrix", "hstack"]


def IM(rows: Iterable[Iterable[int]]) -> IMatrix:
    """Build an immutable matrix from rows"""
    return IMatrix(Matrix(rows).as_immutable())


def IV(components: Iterable[int]) -> IMatrix:
    """Build an immutable vector (immutable matrix with one column)"""
    return IM([[x] for x in components])


def DM(rows: list[list], domain) -> 'DMatrix':
    return DMatrix.from_list(rows, domain)


def DV(components: list, domain) -> 'DMatrix':
    """Build a domain vector (domain matrix with one column)"""
    return DM([[x] for x in components], domain)


def hstack(*args) -> Matrix:
    mats = []
    for m in args:
        if isinstance(m, _DMatrix):
            m = m.to_Matrix()
        mats.append(m)

    if len(mats) == 0:
        return Matrix([])

    n_rows = mats[0].rows
    total_cols = 0
    for m in mats:
        assert m.rows == n_rows
        total_cols += m.cols

    if n_rows == 0 or total_cols == 0:
        return Matrix.zeros(n_rows, total_cols)

    rows = []
    for i in range(n_rows):
        new_row = []
        for m in mats:
            new_row.extend(m.row(i))
        rows.append(new_row)
    return Matrix(rows)


class DMatrix(_DMatrix):
    @classmethod
    def from_list(cls, rows, domain):
        nrows = len(rows)
        ncols = 0 if not nrows else len(rows[0])
        domain_rows = [[domain(e) for e in row] for row in rows]
        return cls(domain_rows, (nrows, ncols), domain)

    def extract_columns(self, indices) -> 'DMatrix':
        return self.extract(list(range(self.shape[0])), indices)

    def columns(self: 'DMatrix') -> list['DMatrix']:
        return [self.extract_columns([i, ]) for i in range(self.shape[1])]

    def __eq__(self, other):
        """The original __eq__ is representation-sensitive. We overwrite it so that it compares values only."""
        return self.to_ddm() == other.to_ddm()

    @classmethod
    def from_Matrix(cls, M: MatrixBase, fmt='sparse', **kwargs):
        assert "domain" in kwargs.keys()
        domain = kwargs["domain"]
        rows = M.tolist()
        domain_rows = [[domain(e) for e in row] for row in rows]
        return cls(domain_rows, (M.rows, M.cols), domain)

    @classmethod
    def static_hstack(cls, A: 'DMatrix', *B: 'DMatrix') -> 'DMatrix':
        assert isinstance(A, DMatrix), type(A)
        res = cls.from_Matrix(hstack(A, *B), domain=A.domain)
        return res

    def __getitem__(self, item):
        from sympy.polys.matrices.domainscalar import DomainScalar
        if isinstance(item, int):
            elem = self.rep.getitem(item // self.shape[1], item % self.shape[1])
            return DomainScalar(self.domain.convert(elem), self.domain)
        else:
            assert len(item) == 2
            return super().__getitem__(item)

    def __str__(self):
        return "NewD" + MatrixBase.__str__(self.to_Matrix())
