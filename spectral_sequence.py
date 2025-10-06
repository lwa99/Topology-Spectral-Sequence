from __future__ import annotations

from element import Bidegree
from page import Page
from utilities import Matrix, Vector, convex_integral_combinations, Poly
from sympy import Symbol
from collections.abc import Iterable


class SpectralSequence:
    def __init__(self, domain, gen: list[Symbol], gen_bidegs, diff_bideg_coef):
        self.gen = gen
        gen_bidegs = Matrix(gen_bidegs)

        assert len(gen) == gen_bidegs.cols
        assert gen_bidegs.rows == 2
        self.generator_bidegs = gen_bidegs

        self.domain = domain  # Base: It's either Z or a field.

        self.pages: list[Page | None] = [None]  # We used None to occupy the 0 index so that page num agrees with index

        # A list of polynomials that we are going to be modded out.
        self.relations: list[Poly] = []

        # A dictionary that maps bidegree to exponents
        self.absolute_bases: dict[Bidegree: tuple[tuple, ...]] = {}

        # The matrix that represents how differential bidegree changes as a function of page number.
        self.diff_bideg_coef = Matrix(diff_bideg_coef)

    @property
    def Z_case(self):
        from sympy import ZZ
        return self.domain == ZZ

    @property
    def gen_num(self):
        return len(self.gen)

    def kill(self, *relations):
        for relation in relations:
            self.relations.append(Poly(relation, *self.gen, domain=self.domain))

    def as_mono(self, exps: Iterable[int]):
        """
        Convert exponent into monomial. For example
        if the generators are {x, y, z}, the tuple (1, 2, 3) is interpreted as x^1y^2z^3
        """
        return Poly.from_dict({tuple(exps): self.domain(1)}, *self.gen, domain=self.domain)

    def get_zero_genset(self, bigrade) -> Matrix:
        """
        Given a bidegree, get a generating set for the zero-space at this bidegree on the first page,
        according to the stored relations.

        Consider all products of relation polynomials with the generators, filter out those having
        the correct bidegree and kill them (by adding them to the basis of the zero-space).

        Note: if the base is Z, this still works but the resulting collection is not a basis. It is a generating
        set of the zero_space.
        """
        if bigrade[1] < 0:
            return Matrix([])

        res = []
        for relation_poly in self.relations:
            abs_bigrade, abs_coordinate = self.get_abs_info(relation_poly)
            skewed_exps = convex_integral_combinations(self.generator_bidegs.row_join(abs_bigrade), bigrade)
            for s_exp in skewed_exps:
                if s_exp[-1] == 0:
                    continue

                temp_bigrade, temp_coordinate = self.get_abs_info(
                    self.as_mono(s_exp[:-1]) * relation_poly ** s_exp[-1]
                )
                print(f"Page 1:\n"
                      f"\tKilled {self.as_mono(s_exp[:-1]) * relation_poly ** s_exp[-1]} in {bigrade}"
                      f"\tas it equals to {self.as_mono(s_exp[:-1])} * "
                      f"{relation_poly} ** {s_exp[-1]} where \t{relation_poly} = 0")
                if temp_coordinate not in res:
                    res.append(temp_coordinate)

        return Matrix.hstack(*res)

    def get_abs_genset(self, bigrade) -> tuple[tuple, ...]:
        """
        Outputs a list of exponents (representing monomials) that gives the correct bidegree.
        """
        if bigrade in self.absolute_bases.keys():
            return self.absolute_bases[bigrade]
        else:
            res = convex_integral_combinations(self.generator_bidegs, bigrade)
            self.absolute_bases[bigrade] = res
            return res

    def get_abs_dimension(self, bigrade: Bidegree):
        """
        In the Z case, this is not a dimension but merely the size of the generating set.
        """
        if bigrade[1] < 0 or (bigrade[1] == 0 and bigrade[0] <= 0):
            return 0
        return len(self.get_abs_genset(bigrade))

    def get_abs_bideg(self, exponent: Iterable[int]) -> Bidegree:
        """
        Compute the bidegree of a monomial (represented by an exponent) by multiplying the
        generator bidegree with exponent.

        Note: generator bidegree are written vertically, so each colum corresponds to a generator.
        """
        return Bidegree(self.generator_bidegs * Matrix(exponent))

    def get_abs_info(self, poly: Poly) -> tuple[Bidegree, Vector]:
        """
        Given a polynomial, compute its bidegree and compute its coordinate on the basis (generating set in Z case).
        """
        abs_bigrade = None
        abs_coordinate = None
        assert not poly.is_zero

        for exponent, coef in poly.terms():
            if abs_bigrade is None:
                abs_bigrade = self.get_abs_bideg(exponent)
            else:
                # Check homogeneity.
                assert abs_bigrade == self.get_abs_bideg(exponent)

            basis = self.get_abs_genset(abs_bigrade)
            cur = []

            # Search through the basis, add a 1 at the corresponding position. Every monomial of the
            # correct degree is guaranteed to occur.
            for n, i in enumerate(basis):
                if i == exponent:
                    cur.append(self.domain(1))
                    cur.extend([self.domain(0)] * (len(basis) - n - 1))
                    break
                else:
                    cur.append(self.domain(0))

            # Multiply the cur (the unit vector) by the coefficient.
            if abs_coordinate is None:
                abs_coordinate = Vector(cur) * coef
            else:
                abs_coordinate += Vector(cur) * coef

        return abs_bigrade, abs_coordinate

    def add_page(self, known_diff: dict = None):
        if known_diff is None:
            known_diff = {}
        page_n = len(self.pages)
        new_page = Page(self, page_n, known_diff, self.diff_bideg_coef * Vector([page_n, 1]))
        self.pages.append(new_page)
        return new_page
