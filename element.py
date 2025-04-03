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


class InHomoPoly:
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
        raise NotImplementedError

    def __sub__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    # def __add__(self, other: 'HomoPoly'):
    #     if len(self.coefMap) == 0:
    #         return deepcopy(other)
    #
    #     # From now on we may assume that self is not the zero polynomial
    #     output = deepcopy(self)
    #     if len(other.coefMap) == 0:
    #         return output
    #     term_deleted = False
    #
    #     for key, scalar in other.coefMap.items():
    #         # If the key is absent, the following line will add the key associated with an empty scalar wrapper.
    #         # Then, the value is replaced by the actual value.
    #         # IF the key is present, the following line will do nothing and return the current value, which will be
    #         # modified accordingly.
    #         # If it becomes 0, the key-value pair will be deleted.
    #         cur_coef = output.coefMap.setdefault(key, self.page.get_scalar(None))
    #         if cur_coef == self.page.get_scalar(None):
    #             cur_coef.update(scalar.val)
    #         else:
    #             cur_coef.increase_by(scalar)
    #
    #             # Now take care of 0 coefficient
    #             if cur_coef.val == 0:
    #                 del output.coefMap[key]
    #                 term_deleted = True
    #
    #     # Now we update bigrades
    #     if output.bigrade is not None:
    #         # If the original bigrade exists but is not equal to that of the other, the bigrade become undefined
    #         if output.bigrade != other.bigrade:
    #             output.bigrade = None
    #     elif term_deleted:
    #         # If the original bigrade did not exist and adding a new polynomial results in deletion of some terms,
    #         # the resulting polynomial may become homogeneous again
    #         output._update_bigrade()
    #
    #     return output
    #
    # def __mul__(self, other):
    #     if len(self.coefMap) == 0:
    #         return deepcopy(self)
    #     if len(other.coefMap) == 0:
    #         return other
    #
    #     # From now on we may assume that both operands are non-zero
    #     output = HomoPoly(self.page)
    #     if self.bigrade is None or other.bigrade is None:
    #         output.bigrade = None
    #     else:
    #         output.bigrade = self.bigrade + other.bigrade
    #
    #     for deg_1, coef_1 in self.coefMap.items():
    #         for deg_2, coef_2 in other.coefMap.items():
    #             output.coefMap[deg_1 + deg_2] = coef_1 * coef_2
    #
    #     return output
    #
    # def _update_bigrade(self):
    #     if len(self.coefMap) == 0:
    #         self.bigrade = Matrix([])
    #         return
    #     _iter = self.coefMap.keys().__iter__()  # get the key iterator
    #     res = self.page.generator_bigrades * _iter.__next__()  # get the first key and calculate bigrade
    #     try:
    #         while True:
    #             if self.page.generator_bigrades * _iter.__next__() != res:
    #                 self.bigrade = None  # Self is non-homogeneous
    #                 return
    #     except StopIteration:
    #         assert isinstance(res, Matrix)
    #         self.bigrade = res
    #         return

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
