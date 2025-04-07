from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page
    from element import HomoElem, Bigrade

from utilities import Matrix, Vector


class Module:
    def __init__(self, page: Page, bigrade: Bigrade, basis: Matrix, ker_basis: Matrix):
        # initialization should only be called by page.getModule
        self.page = page

        """
        The basis and ker_basis here are represented in the standard basis associated with the bigrade.
        """
        self.page = page
        self.bigrade = bigrade

        # Calculate the actual basis.
        # combined: Matrix = Matrix.hstack(ker_basis, basis, Matrix.eye(basis.rows))
        # rref, pivots = combined.rref()
        #
        # ker_basis_idx = []
        # sp_basis_idx = []
        # for i in pivots:
        #     if i < ker_basis.cols:
        #         ker_basis_idx.append(i)
        #     elif i < ker_basis.cols + basis.cols:
        #         sp_basis_idx.append(i - ker_basis.cols)
        ker_basis_idx, sp_basis_idx, self.basis_inv = Matrix.multi_reduction(ker_basis, basis)

        self.sp_basis: Matrix = basis[:, sp_basis_idx]
        self.ker_basis: Matrix = ker_basis[:, ker_basis_idx]

        if self.sp_basis.cols + self.ker_basis.cols > basis.cols:
            raise ValueError("The kernel is not in the span")

    @property
    def basis(self):
        return self.ker_basis.row_join(self.sp_basis)

    @property
    def dim(self):
        return self.page.ss.get_abs_dimension(self.bigrade)

    def __contains__(self, e: HomoElem):
        if e.page != self.page:
            return False
        if e.isZero():
            return True
        return e.bigrade == self.bigrade

    def classify(self, vec: Vector):
        if self.dim == 0:
            if vec.is_zero_matrix:
                return 0
            else:
                return 2
        assert len(vec) == self.dim

        indicator = self.basis_inv * vec
        i = len(indicator) - 1
        while indicator[i] == self.page.ss(0) and i >= 0:
            i -= 1

        if i >= self.ker_basis.cols + self.sp_basis.cols:
            return 2  # not spanned by the basis
        elif i >= self.ker_basis.cols:
            return 1  # spanned by the basis but not the kernel part
        else:
            return 0  # spanned by the kernel basis


if __name__ == "__main__":
    from sympy import pprint
    _basis = Matrix([
        [1,  2],
        [-1, 3],
        [0,  0]
    ])
    _ker_basis = Matrix([
        [1],
        [1],
        [0]
    ])
    m = Module(None, None, _basis, _ker_basis)
    pprint(m.ker_basis)
    pprint(m.sp_basis)
    pprint(m.basis_inv)
    print(m.classify(Vector([3, 3, 1])))
