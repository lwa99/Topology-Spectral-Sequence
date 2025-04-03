"""
There are two types of elements in a page: inhomogeneous polynomial and homogeneous polynomial.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page

from utilities import Vector, Polynomial


class Bigrade(Vector):
    def __repr__(self):
        return f"Bigrade {self.tolist()}"


class InHomoElem:
    pass


class HomoElem:
    """
    Homogeneous polynomials are used to represent homogeneous elements in a page.

    Absolute bigrade: bigrade calculated from the exponent
    Bigrade: the actual bigrade after taking account of the kernel

    Part 1: Initializing



    Zero polynomial:
    Bigrade: None
    Coordinate: None
    """

    def __init__(self,
                 page: Page,
                 poly: Polynomial = None,
                 abs_bigrade: Vector = None,
                 abs_coordinate: Vector = None):
        """
        Two Modes:

        1. Coefficient Map Mode.
        2. Absolute Coordinate Mode
        """
        self.page = page
        ss = page.ss

        if poly is None:
            assert abs_bigrade is not None
            assert abs_coordinate is not None

            # build the coef_map from coordinate
            poly = Polynomial()
            actual_basis = ss.get_abs_basis(abs_bigrade)
            for i, exponent in enumerate(actual_basis):
                poly[exponent] = abs_coordinate[i]

        else:
            assert abs_bigrade is None
            assert abs_coordinate is None

            if len(poly) == 0:
                self.bigrade: Bigrade | None = None
                self.coordinate: Vector | None = None
                self.poly: Polynomial = Polynomial()

            # Get absolute bigrade and coordinate
            abs_bigrade, abs_coordinate = ss.get_abs_info(poly)
        # Step 2: Determine if it is in the kernel and set actual bigrade.
        subspace = page.get_subspace(abs_bigrade)
        r = subspace.classify(abs_coordinate)
        if r == 2:
            raise ValueError("Element does not exist")
        elif r == 1:
            self.bigrade = abs_bigrade
            self.coordinate = abs_coordinate
            self.poly = poly
        else:
            self.bigrade = None
            self.coordinate = None
            self.poly = Polynomial()

    def isZero(self):
        return self.bigrade is None

    def __add__(self, other):
        return HomoElem(self.page, self.poly + other.poly)

    def __sub__(self, other):
        return HomoElem(self.page, self.poly - other.poly)

    def __mul__(self, other):
        return HomoElem(self.page, self.poly * other.poly)


    def __str__(self):
        output = ""
        if len(self.poly) == 0:
            return "zero polynomial (element)"
        for key, value in self.poly.items():
            output += str(value.val)
            for i, degree in enumerate(key):
                output += f"({self.page.ss.generators[i]}^{degree})"
            output += " + "
        output = output[:-2] + f" (mod {self.page.ss.c})"
        return output


if __name__ == "__main__":
    pass
