from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page
    from element import HomoElem

from utilities import Matrix, Vector


class Module:
    def __init__(self, page: Page, bigrade: Vector, basis: Matrix, ker_basis: Matrix):
        # initialization should only be called by page.getModule
        self.page = page

        """
        The basis and ker_basis here are represented in the standard basis associated with the bigrade.
        """
        self.page = page
        self.bigrade = bigrade

        # Calculate the actual basis.
        combined: Matrix = Matrix.hstack(ker_basis, basis, Matrix.eye(basis.rows))
        rref, pivots = combined.rref()

        ker_basis_idx = []
        sp_basis_idx = []
        for i in pivots:
            if i < ker_basis.cols:
                ker_basis_idx.append(i)
            elif i < ker_basis.cols + basis.cols:
                sp_basis_idx.append(i - ker_basis.cols)

        self.sp_basis: Matrix = basis[:, sp_basis_idx]
        self.ker_basis: Matrix = ker_basis[:, ker_basis_idx]

        if self.sp_basis.cols + self.ker_basis.cols > basis.cols:
            raise ValueError("The kernel is not in the span")

        self.basis_inv: Matrix = rref[:, (ker_basis.cols + basis.cols):]

    def __contains__(self, e: HomoElem):
        if e.page != self.page:
            return False
        if e.isZero():
            return True
        return e.bigrade == self.bigrade

    def classify(self, vec: Vector):
        assert len(vec) == self.basis_inv.cols

        indicator = self.basis_inv * vec
        i = len(indicator) - 1
        while indicator[i] == self.page.ss(0):
            i -= 1

        if i == len(indicator) - 1 and i >= self.ker_basis.cols + self.sp_basis.cols:
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
