from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page

from utilities import Vector, Matrix
from sympy import Poly


class Bidegree(Vector):
    def __repr__(self):
        return f"Bigrade {self.tolist()}"


class HomoElem:
    """
    Homogeneous polynomials are used to represent homogeneous elements in a page.

    Absolute bidegree: bidegree calculated from the exponent
    Bidegree: the actual bidegree after taking account of the kernel

    Part 1: Initializing

    Zero polynomial:
    Bigrade: None
    Coordinate: None
    """

    def __init__(self,
                 page: Page,
                 expr = None,
                 abs_bideg: Vector = None,
                 abs_coordinate: Vector = None):
        """
        Two Modes:

        1. Polynomial Mode.
        2. Absolute Coordinate Mode
        """
        self.page = page
        ss = page.ss

        # Step 1: Get abs_info
        if abs_bideg is None:
            assert abs_coordinate is None
            assert expr is not None
            poly = Poly(expr, *ss.gen, domain=ss.domain)
            if poly.total_degree() <= 0:
                self.bidegree: Bidegree | None = None
                self.coordinate: Vector | None = None
                self.poly = poly
                return
            abs_bideg, abs_coordinate = ss.get_abs_info(poly)

        # Step 2: Determine if it is in the kernel and set actual bidegree.
        subspace = page.get_module(abs_bideg)
        r = subspace.classify(abs_coordinate)
        if r == 2:
            raise ValueError(f"Element {expr.__repr__()} "
                             f"(Bigrade: {abs_bideg}) does not exist in page {self.page.page_num}.")
        elif r == 1:
            self.bidegree = abs_bideg
            self.coordinate = abs_coordinate
        else:
            self.bidegree: Bidegree | None = None
            self.coordinate: Vector | None = None
            self.poly = Poly(0, *ss.gen, domain=ss.domain)
            return

        # Step 3: Reconstruct Polynomial in the case that we have a non-trivial element.
        if expr is None:
            assert abs_bideg is not None
            assert abs_coordinate is not None

            # build the polynomial from coordinate
            poly_dict = {}
            actual_basis = ss.get_abs_basis(abs_bideg)
            for i, exponent in enumerate(actual_basis):
                poly_dict[exponent] = abs_coordinate[i]
            self.poly = Poly.from_dict(poly_dict, *ss.gen, domain=ss.domain)
        else:
            self.poly = Poly(expr, *ss.gen, domain=ss.domain)

    def isZero(self):
        return self.bidegree is None

    def __add__(self, other):
        return HomoElem(self.page, self.poly + other.poly)

    def __sub__(self, other):
        return HomoElem(self.page, self.poly - other.poly)

    def __mul__(self, other):
        return HomoElem(self.page, self.poly * other.poly)

    def __eq__(self, other):
        return (self - other).isZero

    def __pow__(self, power, modulo=None):
        if power == 0:
            return HomoElem(self.page, expr="1")
        res = HomoElem(self.page, self.poly)
        for _ in range(power - 1):
            res *= HomoElem(self.page, self.poly)
        return res

    def __hash__(self):
        if self.isZero():
            return 0
        return Matrix(self.coordinate.col_join(self.bidegree)).__hash__()

    def divides(self, other):
        """
        return: whether self divides other, considering the relations.
        """

    def __repr__(self):
        return str(self.poly.as_expr())


if __name__ == "__main__":
    pass
