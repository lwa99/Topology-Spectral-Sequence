from utilities import Matrix
from page import Page
from element import HomoElem


class Differential:
    def __init__(self, page: Page, dx, dy, dz):
        # TODO: allow the differential to be defined via different basis.

        self.page = page
        self.diff_matrix = [dx, dy, dz]
        assert len(self.diff_matrix) == self.page.gen_num  # Ensure diff_matrix is valid

    def calculate(self, exponent_matrix: Matrix) -> HomoElem:
        """
        Compute the differential of a monomial.

        :param exponent_matrix: matrix representing the monomial exponents.
        :return: differential of the monomial
        """
        assert exponent_matrix.shape == (self.page.gen_num, 1)  # Ensure correct shape

        differential_result = HomoElem(self.page)

        for i in range(self.page.gen_num):
            exponent = exponent_matrix[i, 0]  # Get the exponent for this generator
            if exponent > 0:
                new_exponent = exponent_matrix.copy()
                new_exponent[i, 0] -= 1  # Reduce power by 1

                # Convert scalar coefficient to Element
                coefficient_scalar = self.page.get_scalar(exponent)  # take the exponent down when differentiating

                differential_term = self.diff_matrix[i]  # get the image of differential at generator
                term = differential_term * HomoElem(
                    self.page,
                    new_exponent,
                    coefficient_scalar
                )

                differential_result += term  # Sum up according to the product rule

        return differential_result

    # def compute_differential_matrix(self, domain_bigrade: Matrix) -> np.ndarray:
    #     """
    #     对于给定的 domain bigrade，利用当前设置好的 differential，
    #     计算 differential 在该向量空间上的矩阵表达。
    #
    #     流程：
    #     1. 枚举出满足 self.generator_bigrades * exponent = domain_bigrade 的所有指数（即 domain basis）。
    #     2. 对于每个基元（以 Element 表示，系数取 1），计算其 differential。
    #     3. 根据第一个非零 differential 的 bigrade 确定 codomain basis。
    #     4. 将 differential 结果展开成 codomain basis 的线性组合，构造出矩阵表示。
    #     """
    #     # 1. 枚举 domain basis
    #     domain_exponents = enumerate_exponents(self, domain_bigrade)
    #     if not domain_exponents:
    #         print("未找到符合 domain bigrade 的基。")
    #         return np.array([])
    #     domain_basis = [Element(self, exp, self.get_scalar(1)) for exp in domain_exponents]
    #
    #     # 2. 对每个 domain basis 元素应用 differential
    #     differential_images = [self.differential.calculate(exp) for exp in domain_exponents]
    #
    #     # 3. 确定 codomain bigrade（取第一个非零 differential image 的 bigrade）
    #     target_bigrade = None
    #     for img in differential_images:
    #         if len(img.coefMap) > 0 and img.bigrade is not None and len(img.bigrade) > 0:
    #             target_bigrade = img.bigrade
    #             break
    #     if target_bigrade is None:
    #         print("所有 differential 结果均为零。")
    #         return np.zeros((0, len(domain_exponents)), dtype=int)
    #
    #     # 4. 枚举 codomain basis
    #     codomain_exponents = enumerate_exponents(self, target_bigrade)
    #     if not codomain_exponents:
    #         print("未找到符合 codomain bigrade 的基。")
    #         return np.array([])
    #
    #     # 构造矩阵：行对应 codomain basis, 列对应 domain basis
    #     M = np.zeros((len(codomain_exponents), len(domain_exponents)), dtype=int)
    #     for j, img in enumerate(differential_images):
    #         for key, scalar in img.coefMap.items():
    #             found = False
    #             for i, exp in enumerate(codomain_exponents):
    #                 if key == exp:
    #                     M[i, j] = scalar.val
    #                     found = True
    #                     break
    #             if not found:
    #                 print(f"警告：differential 中的项 {key} 不在 codomain basis 内。")
    #     return M


if __name__ == "__main__":
    # 创建 Page 实例
    generators = ["x", "y", "z"]
    generator_bigrades = Matrix([[1, 0, 3], [4, -3, 6]])  # 生成元的双次数
    c = 7  # 有限域特征数

    _page = Page(generators, generator_bigrades, c, 1)

    # 生成 d(x), d(y), d(z) 作为 Element
    _dx = HomoElem(_page, Matrix([[2], [0], [0]]), Scalar(c, 1))
    _dy = HomoElem(_page, Matrix([[0], [2], [0]]), Scalar(c, 1))
    _dz = HomoElem(_page, Matrix([[0], [0], [2]]), Scalar(c, 1))  # d(z) 对应于 x^0 y^0 z^2

    # 创建 Differential 实例
    differential = Differential(_page, _dx, _dy, _dz)

    # 选取 exponent_matrix，表示 x^1 y^2 z^3
    _exponent_matrix = Matrix([[1], [2], [3]])

    # 计算微分
    result = differential.calculate(_exponent_matrix)

    # 打印结果
    print("微分计算结果: ", result)
