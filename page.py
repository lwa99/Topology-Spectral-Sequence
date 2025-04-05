from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from spectral_sequence import SpectralSequence

from utilities import Matrix, Polynomial, degree_generator
from sortedcontainers import SortedDict
from differential import Differential
from element import HomoElem
from sympy import div
from module import Module


class Page:
    def __init__(self, ss, page_num, d: Differential):
        self.ss: SpectralSequence = ss
        self.page_num = page_num
        self.subspaces = SortedDict([])
        self.d: Differential = d

    def get_module(self, bigrade):
        if bigrade in self.subspaces:
            return self.subspaces[bigrade]
        else:
            output = self.generate_module(bigrade)
            self.subspaces[bigrade] = output
            return output

    def generate_module(self, bigrade) -> Module:
        # This is the function that handles generates subspaces of new pages.
        if self.page_num == 1:
            return Module(self, bigrade,
                          Matrix.eye(self.ss.get_abs_dimension(bigrade)),
                          self.ss.get_ker_basis(bigrade)
                          )

        """
        计算 A2 的核（kernel）以及 A1 的像（image）。

        参数：
        - M1: 模块 M1，具有一个已知的 bigrade。
        - differential_bigrade: 微分的 bigrade，用于计算 M0。
        - get_matrix: 一个方法，输入 bigrade，返回相应的矩阵。

        返回：
        - kernel_A2: A2 的核（null space）。
        - image_A1: A1 的像（column space）。
        """
        # 计算 prev 的 bigrade
        prev_bigrade = bigrade - self.d.d_bigrade  # 直接用向量减法

        # 获取 A1 和 A2
        matrix_prev = self.d.get_matrix(prev_bigrade)  # A1: M0 -> M1 的映射矩阵
        matrix_next = self.d.get_matrix(bigrade)  # A2: M1 -> M_next 的映射矩阵

        # 计算 A2 的核（解 Ax=0）
        kernel_next = matrix_next.nullspace()

        # 计算 A1 的像（列空间）
        image_prev = matrix_prev.columnspace()

        return Module(self, bigrade, image_prev, kernel_next)

    def find_kernels_for_division(self,
                                  a: Polynomial,
                                  c: Polynomial,
                                  bigrade):
        # 大体method没有问题，需要确认有效性，和mod prime_char
        """
        找到所有 (k1, k2)，使得 (c + k1) / (a + k2) 整除成立。

        参数：
        - a, c: SymPy 多项式
        - ker_basis_matrix: kernel 的 basis（矩阵）
        - page: Page 对象，用于构造 HomoElem
        - bigrade: 当前处理的 bigrade
        - prime_char: 有限域的特征（默认是 3）

        返回：
        - solutions: 满足整除条件的三元组 (k1, k2, b)
        """
        kernels = []
        module = self.get_module(bigrade)
        ker_basis = module.ker_basis
        char = self.ss.c
        for degree_vector in degree_generator(ker_basis, char):
            try:
                elem = HomoElem(self, poly=None, abs_bigrade=bigrade, abs_coordinate=degree_vector)
                if not elem.isZero():
                    kernels.append(elem.poly)
            except ValueError:
                continue

        solutions = []
        for k1 in kernels:
            for k2 in kernels:
                numerator = c + k1
                denominator = a + k2
                q, r = div(numerator, denominator, domain=self.ss.ff)
                if r == 0:
                    solutions.append((k1, k2, q))

        return solutions

    # def test_zero_homo_poly(self, coef_map: SortedDict):
    #     assert len(Vector) == len(self.ss.generators)
    #     first_exponent = coef_map.keys().__iter__().__next__()
    #     bigrade = self.ss.get_bigrade(first_exponent)
    #     subspace = self.get_subspace(bigrade)
    #     return subspace.kernelContains(self.ss.get_std_coordinate(coef_map))


if __name__ == "__main__":
    pass
