from __future__ import annotations

from element import Bigrade
from page import Page
from utilities import Matrix, Vector, convex_integral_combinations
from sympy import Poly, Symbol, Expr
from collections.abc import Iterable


class SpectralSequence:
    def __init__(self, gen: list[Symbol], generator_bigrades, domain):
        self.gen: list[Poly] = [Poly(g, *gen, domain=domain) for g in gen]
        generator_bigrades = Matrix(generator_bigrades)

        assert len(gen) == generator_bigrades.cols
        assert generator_bigrades.rows == 2
        self.generator_bigrades = generator_bigrades

        self.domain = domain  # Base Field

        self.pages: list[Page | None] = [None]  # We used None to occupy the 0 index so that page num agrees with index
        self.relations: list[Poly] = []

        # A dictionary that maps bigrades to exponents
        self.absolute_bases: dict[Bigrade: tuple[tuple, ...]] = {}

    @property
    def gen_num(self):
        return len(self.gen)

    # def kill(self, *polys: Polynomial):
    #     self.relations.extend(polys)

    def add_relation(self, relation: Expr):
        self.relations.append(Poly(relation, *self.gen, domain=self.domain))

    def as_mono(self, exps: Iterable[int]):
        """
        Convert exponent into monomial
        """
        return Poly.from_dict({tuple(exps): self.domain(1)}, *self.gen, domain=self.domain)

    def get_ker_basis(self, bigrade) -> Matrix:
        """
        Calculate the first page kernel basis at a given bigrade from stored relations
        """
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
                print(f"Page 1:"
                      f"Killed {self.as_mono(s_exp[:-1]) * relation_poly ** s_exp[-1]} in {bigrade} as it equals to"
                      f"{self.as_mono(s_exp[:-1])} * {relation_poly} ** {s_exp[-1]} where {relation_poly} is killed "
                      f"initially.")
                if temp_coordinate not in res:
                    res.append(temp_coordinate)
        return Matrix.hstack(*res)

    def get_abs_basis(self, bigrade) -> tuple[tuple, ...]:
        if bigrade in self.absolute_bases.keys():
            return self.absolute_bases[bigrade]
        else:
            res = convex_integral_combinations(self.generator_bigrades, bigrade)
            self.absolute_bases[bigrade] = res
            return res

    def get_abs_dimension(self, bigrade: Bigrade):
        if bigrade[1] < 0 or (bigrade[1] == 0 and bigrade[0] <= 0):
            return 0
        return len(self.get_abs_basis(bigrade))

    def get_abs_bigrade(self, exponent: Iterable[int]) -> Bigrade:
        return Bigrade(self.generator_bigrades * Matrix(exponent))

    def get_abs_info(self, poly: Poly) -> tuple[Bigrade, Vector]:
        abs_bigrade = None
        abs_coordinate = None
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
                abs_coordinate = Vector(cur) * coef
            else:
                abs_coordinate += Vector(cur) * coef

        return abs_bigrade, abs_coordinate

    # def get_std_coordinate(self, coef_map: SortedDict) -> Vector | None:
    #     bigrade = None
    #     res = None
    #     assert len(coef_map) > 0
    #     for exponent, coef in coef_map.items():
    #         if bigrade is None:
    #             bigrade = self.get_abs_bigrade(exponent)
    #         else:
    #             assert bigrade == self.get_abs_bigrade(exponent)
    #
    #         basis = self.get_actual_basis(bigrade)
    #         cur = []
    #         for n, i in enumerate(basis):
    #             if i == exponent:
    #                 cur.append(1)
    #                 cur.extend([0] * (len(basis) - n - 1))
    #                 break
    #             else:
    #                 cur.append(0)
    #         if res is None:
    #             res = Vector(cur) * coef
    #         else:
    #             res += Vector(cur) * coef
    #     return res
