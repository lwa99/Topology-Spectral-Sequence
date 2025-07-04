from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from spectral_sequence import SpectralSequence

from utilities import Matrix, degree_generator
from sortedcontainers import SortedDict
from differential import Differential
from element import Bidegree
from sympy import div
from module import Module


class Page:
    def __init__(self, ss, page_num, io_pairs: dict, d_bigrade):
        self.ss: SpectralSequence = ss
        self.page_num = page_num
        self.subspaces = SortedDict([])
        self.d = Differential(self, io_pairs, Bidegree(d_bigrade))

    def get_module(self, bigrade: Bidegree):
        if bigrade in self.subspaces:
            return self.subspaces[bigrade]
        else:
            output = self.generate_module(bigrade)
            self.subspaces[bigrade] = output
            return output

    def generate_module(self, bigrade) -> Module:
        if self.page_num == 1:
            return Module(self, bigrade,
                          Matrix.eye(self.ss.get_abs_dimension(bigrade)),
                          self.ss.get_ker_basis(bigrade)
                          )

        prev_page = self.ss.pages[self.page_num - 1]
        # 计算 prev 的 bidegree
        prev_bigrade = bigrade - prev_page.d.d_bidegree  # 直接用向量减法

        # 获取 A1 和 A2
        matrix_prev = prev_page.d.get_matrix(prev_bigrade, debug=self.page_num == 4)  # A1: M0 -> M1 的映射矩阵
        matrix_next = prev_page.d.get_matrix(bigrade)  # A2: M1 -> M_next 的映射矩阵

        # 计算 A2 的核（解 Ax=0）
        kernel_next = matrix_next.nullspace()

        # 计算 A1 的像（列空间）
        image_prev = matrix_prev.columnspace()

        return Module(self, bigrade,
                      Matrix.hstack(*kernel_next),
                      Matrix.hstack(*image_prev, prev_page[bigrade].ker_basis))

    def __getitem__(self, item):
        return self.get_module(Bidegree(item))

    # def find_kernels_for_division(self,
    #                               a: Polynomial,
    #                               c: Polynomial,
    #                               bidegree):
    #     # 大体method没有问题，需要确认有效性，和mod prime_char
    #     """
    #     找到所有 (k1, k2)，使得 (c + k1) / (a + k2) 整除成立。
    #
    #     参数：
    #     - a, c: SymPy 多项式
    #     - ker_basis_matrix: kernel 的 basis（矩阵）
    #     - page: Page 对象，用于构造 HomoElem
    #     - bidegree: 当前处理的 bidegree
    #     - prime_char: 有限域的特征（默认是 3）
    #
    #     返回：
    #     - solutions: 满足整除条件的三元组 (k1, k2, b)
    #     """
    #     kernels = []
    #     module = self.get_module(bidegree)
    #     ker_basis = module.ker_basis
    #     char = self.ss.c
    #     for degree_vector in degree_generator(ker_basis, char):
    #         try:
    #             elem = HomoElem(self, poly=None, abs_bigrade=bidegree, abs_coordinate=degree_vector)
    #             if not elem.isZero():
    #                 kernels.append(elem.poly)
    #         except ValueError:
    #             continue
    #
    #     solutions = []
    #     for k1 in kernels:
    #         for k2 in kernels:
    #             numerator = c + k1
    #             denominator = a + k2
    #             q, r = div(numerator, denominator, domain=self.ss.bf)
    #             if r == 0:
    #                 solutions.append((k1, k2, q))
    #
    #     return solutions

    # def test_zero_homo_poly(self, coef_map: SortedDict):
    #     assert len(Vector) == len(self.ss.generators)
    #     first_exponent = coef_map.keys().__iter__().__next__()
    #     bidegree = self.ss.get_bigrade(first_exponent)
    #     subspace = self.get_subspace(bidegree)
    #     return subspace.kernelContains(self.ss.get_std_coordinate(coef_map))


if __name__ == "__main__":
    pass
