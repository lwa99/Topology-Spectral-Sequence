from tools import DOArray
from element import Element
from page import Page
import numpy as np


class Initializer:
    def __init__(self, page, bigrade):
        self.page = page
        self.bigrade = np.array(bigrade)
        self.basis = []
        self._generate_basis([], 0)

    def _generate_basis(self, current_exp, idx):
        """
        递归生成所有可能的指数组合，使得 generator_bigrades @ exponents == bigrade
        """
        if idx == self.page.gen_num:  # 所有变量的指数都已经选定
            current_exp = np.array(current_exp)
            if np.all(current_exp @ self.page.generator_bigrades == self.bigrade):
                self.basis.append(DOArray(current_exp))
            return

        # 计算当前变量的指数上限
        max_exp = self.bigrade // self.page.generator_bigrades[:, idx]
        max_exp = np.min(max_exp)  # 取所有限制中的最小值，防止指数超限

        # 递归枚举当前变量的指数
        for exp in range(max(0, max_exp) + 1):
            self._generate_basis(current_exp + [exp], idx + 1)

    def get_basis(self):
        return self.basis


if __name__ == "__main__":
    # 生成测试 Page
    generators = ["x", "y"]
    generator_bigrades = DOArray([[1, 2], [3, 4]])  # 2x2 矩阵，每列是一个 generator 的 bigrade
    c = 7  # 质数特征

    _page = Page(generators, generator_bigrades, c)

    # 设定要测试的 bigrade
    _bigrade = [7, 10]

    # 初始化并计算 basis
    initializer = Initializer(_page, _bigrade)
    basis = initializer.get_basis()

    # 打印 basis 结果
    print("Computed Basis:")
    for monomial in basis:
        print(monomial)
