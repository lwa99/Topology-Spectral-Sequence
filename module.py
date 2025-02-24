from tools import Matrix, pprint
from page import Page
import numpy as np


class Initializer:
    def __init__(self, page: Page, bigrade: Matrix):
        self.page = page
        assert bigrade.shape == (2, 1)
        self.bigrade = bigrade
        self.basis = []
        self._calculate_basis()

    def _generate_basis(self, current_exp, idx):
        """
        递归生成所有可能的指数组合，使得 generator_bigrades @ exponents == bigrade
        """
        if idx == self.page.gen_num:  # 所有变量的指数都已经选定
            current_exp = np.array(current_exp)
            if np.all(current_exp @ np.array(self.page.generator_bigrades) == self.bigrade):
                self.basis.append(Matrix(current_exp))
            return

        # 计算当前变量的指数上限
        max_exp = self.bigrade // np.array(self.page.generator_bigrades)[:, idx]
        max_exp = np.min(max_exp)  # 取所有限制中的最小值，防止指数超限

        # 递归枚举当前变量的指数
        for exp in range(max(0, max_exp) + 1):
            self._generate_basis(current_exp + [exp], idx + 1)

    def _calculate_basis(self):
        bg = self.page.generator_bigrades

        # Classical Adjugate of a 2x2 matrix
        # inv = Matrix([[bg[1, 1], -bg[0, 1]], [-bg[1, 0], bg[0, 0]]])
        # b = inv * self.bigrade

        # Now, we want to find the bound of the variables
        bounds = [-1] * self.page.gen_num
        skipped_index = []
        for j in range(self.page.gen_num):
            print(j)
            if bg[0, j] > 0:
                bounds[j] = self.bigrade[0] // bg[0, j]
                if bg[1, j] > 0 and self.bigrade[1] // bg[1, j] < bounds[j]:
                    bounds[j] = self.bigrade[1] // bg[1, j]
                continue
            skipped_index.append(j)

        # Now, we are left with those generators with 0 in the first grade
        if len(skipped_index) > 0:
            cap = self.bigrade[1]
            for j in range(self.page.gen_num):
                if j in skipped_index:
                    continue
                if bg[1, j] * bg[1, skipped_index[0]] < 0:
                    cap -= bg[1, j] * bounds[j]

            for j in skipped_index:
                if bg[1, j] * bg[1, skipped_index[0]] <= 0:
                    raise RuntimeError
                bounds[j] = cap // bg[1, j]
                if bounds[j] < 0:
                    return  # there is no possible exponent that can give the bigrade.

        print(bounds)

    def get_basis(self):
        return self.basis


if __name__ == "__main__":
    # 生成测试 Page
    generators = ["x", "y", "z"]
    generator_bigrades = Matrix([[1, 0, 3], [4, -3, 6]])  # 2x2 矩阵，每列是一个 generator 的 bigrade
    c = 7  # 质数特征
    print("Generator bigrades:")
    pprint(generator_bigrades)

    _page = Page(generators, generator_bigrades, c)

    # 设定要测试的 bigrade
    _bigrade = Matrix([7, 10])
    print("Target:")
    pprint(_bigrade)

    # 初始化并计算 basis
    initializer = Initializer(_page, _bigrade)
    basis = initializer.get_basis()

    # 打印 basis 结果
    print("Computed Basis:")
    for monomial in basis:
        print(monomial)
