from __future__ import annotations

from .element import Bidegree
from page_and_module import Page
from utilities import convex_integral_combinations, Poly
from matrices import *
from sympy import Symbol
from collections.abc import Iterable


class SpectralSequence:
    def __init__(self, domain, gen: list[Symbol], generator_bideg: list[list], diff_bideg_coef):
        self.gen = gen
        generator_bideg = IM(generator_bideg)

        assert len(gen) == generator_bideg.cols  # each generator corresponds to one column
        assert generator_bideg.rows == 2
        self.generator_bigrades = generator_bideg

        self.domain = domain  # Base Field

        self.pages: list[Page | None] = [None]  # We used None to occupy the 0 index so that page num agrees with index
        self.relations: list[Poly] = []

        # A dictionary that maps bidegree to exponents
        self.absolute_bases: dict[Bidegree: tuple[tuple, ...]] = {}

        self.diff_bideg_coef = IM(diff_bideg_coef)

    @property
    def gen_num(self):
        return len(self.gen)

    def kill(self, *relations):
        for relation in relations:
            self.relations.append(Poly(relation, *self.gen, domain=self.domain))

    def as_mono(self, exps: Iterable[int]):
        """
        Convert exponent into monomial
        """
        return Poly.from_dict({tuple(exps): self.domain(1)}, *self.gen, domain=self.domain)

    def get_ker_basis(self, bigrade) -> list[DMatrix]:
        """
        Calculate the first page kernel basis at a given bidegree from stored relations
        """
        if bigrade[1] < 0:
            return DM([], self.domain).columns()

        res = []
        for relation_poly in self.relations:
            abs_bigrade, abs_coordinate = self.get_abs_info(relation_poly)
            skewed_exps = convex_integral_combinations(self.generator_bigrades.row_join(abs_bigrade), bigrade)
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
        return res

    def get_abs_basis(self, bigrade) -> tuple[tuple, ...]:
        if bigrade in self.absolute_bases.keys():
            return self.absolute_bases[bigrade]
        else:
            res = convex_integral_combinations(self.generator_bigrades, bigrade)
            self.absolute_bases[bigrade] = res
            return res

    def get_abs_dimension(self, bigrade: Bidegree):
        if bigrade[1] < 0 or (bigrade[1] == 0 and bigrade[0] < 0):
            return 0
        if bigrade[1] == bigrade[0] == 0:
            return 1  # the scalar space
        return len(self.get_abs_basis(bigrade))

    def get_abs_bigrade(self, exponent: Iterable[int]) -> Bidegree:
        return Bidegree(self.generator_bigrades * IV(exponent))

    def get_abs_info(self, poly: Poly) -> tuple[Bidegree | None, DMatrix | None]:
        abs_bigrade = None
        abs_coordinate = None
        if poly.total_degree() == 0:
            assert len(poly.terms()) == 1
            if poly.terms()[0][1] == self.domain.zero:
                return None, None
            return IV([0, 0]), DV([poly.terms()[0][1]], self.domain)
        assert not poly.is_zero
        for exponent, coef in poly.terms():
            if abs_bigrade is None:
                abs_bigrade = self.get_abs_bigrade(exponent)
            else:
                assert abs_bigrade == self.get_abs_bigrade(exponent)

            basis = self.get_abs_basis(abs_bigrade)
            cur = []
            for n, i in enumerate(basis):
                if i == exponent:
                    cur.append(self.domain(1))
                    cur.extend([self.domain(0)] * (len(basis) - n - 1))
                    break
                else:
                    cur.append(self.domain(0))
            if abs_coordinate is None:
                abs_coordinate = DV(cur, self.domain) * coef
            else:
                abs_coordinate += DV(cur, self.domain) * coef

        return abs_bigrade, abs_coordinate

    def add_page(self, known_diff: dict = None):
        if known_diff is None:
            known_diff = {}
        page_n = len(self.pages)
        new_page = Page(self, page_n, known_diff, self.diff_bideg_coef * IV([page_n, 1]))
        self.pages.append(new_page)
        return new_page


if __name__ == "__main__":
    from sympy.abc import symbols
    from sympy import GF, ZZ
    from element import HomoElem, HomoCollection
    a, t = symbols("a t")
    ss = SpectralSequence(
        ZZ,
        [a, t],
        [[3, 3], [0, 0]],
        [[1, 0],
         [-1, 1]]
    )

    # print("abs_dim", ss.get_abs_dimension(IV([3, 0])))
    # ss.kill(a-t)
    # p1 = ss.add_page({a: 0, t: 0})
    # target = HomoElem(p1, 2*t)
    #
    # m = HomoElem(p1, a)
    # r = p1.divide(m, target)
    print("abs_basis", ss.get_abs_basis(IV([6, 0])))
    p1 = ss.add_page({a: 1, t: 1})
    a_2 = HomoElem(p1, a**2)
    at = HomoElem(p1, a*t)
    t_2 = HomoElem(p1, t**2)

    c_1 = a_2 + t_2
    c_2 = a_2 + at
    c_3 = a_2 - t_2

    r_1 = HomoElem(p1, 2*a + 2*t)
    r_2 = HomoElem(p1, 3*a + t)
    r_3 = HomoElem(p1, 2*a - 2*t)

    I = HomoCollection([c_1, c_2, c_3])
    d_I = HomoCollection([r_1, r_2, r_3])
    span = p1[IV([6, 0])].get_diff_span(I, d_I)
    print(span.elems)

