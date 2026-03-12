from __future__ import annotations
from typing import TYPE_CHECKING

from utilities import verify
from warnings import warn
from matrices import *
from sympy import Poly, Expr
from collections.abc import Iterable
from snf import SNFMatrix

if TYPE_CHECKING:
    from page_and_module import Page

Bidegree = IMatrix


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
                 expr=None,
                 *,
                 abs_bideg: IMatrix = None,
                 abs_coordinate: DMatrix = None):
        """
        Two Modes:

        1. Polynomial Mode.
        2. Absolute Coordinate Mode
        """
        self.page = page
        ss = page.ss

        if expr is None:
            assert abs_bideg is not None
            assert abs_coordinate is not None

            self.bidegree = abs_bideg
            self.coordinate = abs_coordinate
            # build the polynomial from coordinate
            poly_dict = {}
            abs_basis = ss.get_abs_basis(abs_bideg)
            for i, exponent in enumerate(abs_basis):
                poly_dict[exponent] = abs_coordinate[i].element
            self.poly = Poly.from_dict(poly_dict, *ss.gen, domain=ss.domain)
        else:
            assert abs_coordinate is None
            assert abs_bideg is None

            self.poly = Poly(expr, *ss.gen, domain=ss.domain)
            self.bidegree, self.coordinate = ss.get_abs_info(self.poly)

    def isZero(self):
        return self.bidegree is None

    def __add__(self, other):
        return HomoElem(self.page, self.poly + other.poly)

    def __sub__(self, other):
        return HomoElem(self.page, self.poly - other.poly)

    def __mul__(self, other):
        if isinstance(other, HomoCollection):
            return other.__rmul__(self)
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
        return hash(tuple(self.coordinate.to_list_flat())) + self.bidegree.__hash__()

    def divides(self, other):
        """
        return: whether self divides other, considering the relations.
        """
        raise NotImplementedError

    def __repr__(self):
        return str(self.poly.as_expr())


class HomoCollection:
    """Collection of HomoElem of the same bidegree"""
    def __init__(self,
                 page: Page | None,
                 bideg: Bidegree | None,
                 elems: Iterable[HomoElem] = None,
                 coords: Iterable[DMatrix] = None):

        if elems is not None:
            self._elems = list(elems)
            self._coords = None

            if len(self._elems) == 0:
                assert page is not None
                assert bideg is not None
                self.bideg = bideg
                self.page = page
            else:
                self.bideg = self._elems[0].bidegree
                self.page = self._elems[0].page
                if bideg is not None:
                    assert bideg == self.bideg
                if page is not None:
                    assert page == self.page
                if verify:
                    for e in self._elems:
                        assert e.page == self.page
                        assert e.bidegree == self.bideg

        else:
            assert page is not None
            assert bideg is not None
            assert coords is not None
            self._coords = list(coords)
            self._elems = None
            self.bideg = bideg
            self.page = page

    @property
    def elems(self) -> list[HomoElem]:
        if self._elems is None:
            self._elems = [HomoElem(self.page, abs_bideg=self.bideg, abs_coordinate=v) for v in self._coords]
        return self._elems

    @property
    def coords(self) -> list[DMatrix]:
        if self._coords is None:
            self._coords = [e.coordinate for e in self._elems]
        return self._coords

    @property
    def is_empty(self):
        return len(self) == 0

    @staticmethod
    def from_elems(elems):
        """Only use when elems is non-empty."""
        assert len(elems) > 0
        return HomoCollection(None, None, elems)

    @staticmethod
    def from_exprs(page, exprs: Iterable[Expr]):
        """Only use when exprs is non-empty."""
        data = [HomoElem(page, expr) for expr in exprs]
        return HomoCollection.from_elems(data)

    @staticmethod
    def from_matrix(page: Page, bideg: Bidegree, M: DMatrix):
        return HomoCollection(page=page, bideg=bideg, coords=M.columns())

    def to_matrix(self) -> DMatrix | None:
        if self.is_empty:
            return None
        return DMatrix.static_hstack(*self.coords)

    def to_SNF_matrix(self) -> SNFMatrix | None:
        if self.is_empty:
            return None
        return SNFMatrix.static_hstack(*self.coords)

    def __len__(self):
        if self._elems is not None:
            return len(self._elems)
        else:
            return len(self._coords)

    def __mul__(self, M: DMatrix):
        S = self.to_matrix()
        if S is None:
            return self
        return HomoCollection.from_matrix(self.page, self.bideg, S * M)

    def __rmul__(self, x: HomoElem):
        assert self.page == x.page
        return HomoCollection(self.page, x.bidegree + self.bideg, elems=[x * y for y in self.elems])

    def join(self, other: HomoCollection):
        assert self.page == other.page and self.bideg == other.bideg
        return HomoCollection(self.page, other.bideg, elems=self.elems + other.elems)


if __name__ == "__main__":
    pass
