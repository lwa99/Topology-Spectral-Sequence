from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from spectral_sequence import SpectralSequence

from utilities import Matrix
from sortedcontainers import SortedDict
from differential import Differential
from element import Bidegree
from module import Module


class Page:
    def __init__(self, ss, page_num, io_pairs: dict, d_bigrade):
        self.ss: SpectralSequence = ss
        self.page_num = page_num
        self.subspaces = SortedDict([])
        self.d = Differential(self, io_pairs, d_bigrade)

    def get_module(self, bigrade: Bidegree):
        if bigrade in self.modules:
            return self.modules[bigrade]
        else:
            output = self.generate_module(bigrade)
            self.modules[bigrade] = output
            return output

    def generate_module(self, bigrade) -> Module:
        # First Page
        if self.page_num == 1:
            return Module(self, bigrade,
                          Matrix.eye(self.ss.get_abs_dimension(bigrade)),
                          self.ss.get_zero_genset(bigrade)
                          )

        prev_page = self.ss.pages[self.page_num - 1]
        # 计算 prev 的 bidegree
        prev_bigrade = bigrade - prev_page.d.d_bidegree  # 直接用向量减法

        # 获取 A1 和 A2
        matrix_prev = prev_page.d.get_matrix(prev_bigrade)  # A1: M0 -> M1 的映射矩阵
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


if __name__ == "__main__":
    pass
