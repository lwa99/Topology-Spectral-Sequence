from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from element import Bigrade
    from page import Page

from utilities import Matrix, Vector, Prime, Polynomial, Exponent, convex_integral_combinations, monomial
from sympy import GF


class SpectralSequence:
    def __init__(self, generators: list[str], generator_bigrades: Matrix, c: int):
        self.generators = generators

        assert len(generators) == generator_bigrades.cols
        assert generator_bigrades.rows == 2
        self.generator_bigrades = generator_bigrades

        assert Prime.is_prime(c)
        self.c = c
        self.ff = GF(c)
        Polynomial.initiate(self.ff, self.generators)

        self.pages: list[Page] = []
        self.relations: list[Polynomial] = []

        # A dictionary that maps bigrades to exponents
        self.actual_bases: dict[Bigrade: tuple[Exponent, ...]] = {}

    def add_relation(self, relation: Polynomial):
        self.relations.append(relation)

    def get_ker_basis(self, bigrade) -> Matrix:
        """
        Calculate the first page kernel basis at a given bigrade from stored relations
        """
        res = []
        for relation in self.relations:
            abs_bigrade, abs_coordinate = self.get_absolute_info(relation)
            skewed_exps = convex_integral_combinations(self.generator_bigrades.row_join(abs_bigrade), bigrade)
            for s_exp in skewed_exps:
                if s_exp[-1] == 0:
                    continue

                temp_bigrade, temp_coordinate = self.get_absolute_info(
                    monomial(s_exp[:-1]) * relation ** s_exp[-1]
                )
                print("Formula:", f"{monomial(s_exp[:-1])} * {relation} ** {s_exp[-1]} = "
                                  f"{monomial(s_exp[:-1]) * relation ** s_exp[-1]}")
                print()
                if temp_coordinate not in res:
                    res.append(temp_coordinate)
        return Matrix.hstack(*res)

    def get_actual_basis(self, bigrade) -> tuple[Exponent, ...]:
        if bigrade in self.actual_bases.keys():
            return self.actual_bases[bigrade]
        else:
            res = tuple([Exponent(v) for v in convex_integral_combinations(self.generator_bigrades, bigrade)])
            self.actual_bases[bigrade] = res
            return res

    def get_abs_dimension(self, bigrade: Bigrade):
        return len(self.get_actual_basis(bigrade))

    def get_abs_bigrade(self, exponent) -> Bigrade:
        return self.generator_bigrades * exponent

    def get_absolute_info(self, coef_map) -> tuple[Bigrade, Vector]:
        print("abs_info called on:", coef_map.__str__())
        abs_bigrade = None
        abs_coordinate = None
        assert len(coef_map) > 0
        for exponent, coef in coef_map.items():
            if abs_bigrade is None:
                abs_bigrade = self.get_abs_bigrade(exponent)
            else:
                assert abs_bigrade == self.get_abs_bigrade(exponent)

            basis = self.get_actual_basis(abs_bigrade)
            cur = []
            for n, i in enumerate(basis):
                if i == exponent:
                    cur.append(self(1))
                    cur.extend([self(0)] * (len(basis) - n - 1))
                    break
                else:
                    cur.append(self(0))
            if abs_coordinate is None:
                abs_coordinate = Vector(cur) * coef
            else:
                abs_coordinate += Vector(cur) * coef

        print("abs_info returns:", abs_bigrade, abs_coordinate)
        return abs_bigrade, abs_coordinate

    def __call__(self, val):
        if isinstance(val, int):
            return self.ff(val)
        elif isinstance(val, dict):
            return Polynomial(dict)
        else:
            return monomial(val)

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
