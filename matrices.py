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
    non_empty_args = []
    for m in args:
        if isinstance(m, _DMatrix):
            m = m.to_Matrix()
        if m.rows != 0 and m.cols != 0:
            non_empty_args.append(m)

    rows = []
    n_rows = non_empty_args[0].rows
    for i in range(n_rows):
        new_row = []
        for m in non_empty_args:
            assert m.rows == n_rows
            new_row.extend(m.row(i))
        rows.append(new_row)
    return Matrix(rows)


class DMatrix(_DMatrix):
    @classmethod
    def _from_original_DMatrix(cls, M: _DMatrix) -> 'DMatrix':
        return DMatrix(M.to_list(), M.shape, M.domain)

    @classmethod
    def from_list(cls, rows, domain):
        res = cls._from_original_DMatrix(super().from_list(rows, domain))
        return res

    def extract_columns(self, indices) -> 'DMatrix':
        return self.extract(list(range(self.shape[0])), indices)

    def columns(self: 'DMatrix') -> list['DMatrix']:
        return [self.extract_columns([i, ]) for i in range(self.shape[1])]

    def __eq__(self, other):
        return self.to_ddm() == other.to_ddm()

    @classmethod
    def from_Matrix(cls, M: MatrixBase, fmt='sparse', **kwargs):
        assert "domain" in kwargs.keys()
        return cls.from_list(M.tolist(), **kwargs)

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
