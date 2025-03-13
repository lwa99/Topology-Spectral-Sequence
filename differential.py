from tools import Matrix
from page import Page
from element import Element
from scalar import Scalar


class Differential:
    def __init__(self, page: Page, dx, dy, dz):
        # TODO: allow the differential to be defined via different basis.

        self.page = page
        self.diff_matrix = [dx, dy, dz]
        assert len(self.diff_matrix) == self.page.gen_num  # Ensure diff_matrix is valid

    def calculate(self, exponent_matrix: Matrix) -> Element:
        """
        Compute the differential of a monomial.

        :param exponent_matrix: matrix representing the monomial exponents.
        :return: differential of the monomial
        """
        assert exponent_matrix.shape == (self.page.gen_num, 1)  # Ensure correct shape

        differential_result = Element(self.page)

        for i in range(self.page.gen_num):
            exponent = exponent_matrix[i, 0]  # Get the exponent for this generator
            if exponent > 0:
                new_exponent = exponent_matrix.copy()
                new_exponent[i, 0] -= 1  # Reduce power by 1

                # Convert scalar coefficient to Element
                coefficient_scalar = self.page.get_scalar(exponent)  # take the exponent down when differentiating

                differential_term = self.diff_matrix[i]  # get the image of differential at generator
                term = differential_term * Element(
                    self.page,
                    new_exponent,
                    coefficient_scalar
                )

                differential_result += term  # Sum up according to the product rule

        return differential_result



if __name__ == "__main__":
    # 创建 Page 实例
    generators = ["x", "y", "z"]
    generator_bigrades = Matrix([[1, 0, 3], [4, -3, 6]])  # 生成元的双次数
    c = 7  # 有限域特征数

    _page = Page(generators, generator_bigrades, c, 1)

    # 生成 d(x), d(y), d(z) 作为 Element
    _dx = Element(_page, Matrix([[2], [0], [0]]), Scalar(c, 1))
    _dy = Element(_page, Matrix([[0], [2], [0]]), Scalar(c, 1))
    _dz = Element(_page, Matrix([[0], [0], [2]]), Scalar(c, 1))  # d(z) 对应于 x^0 y^0 z^2

    # 创建 Differential 实例
    differential = Differential(_page, _dx, _dy, _dz)

    # 选取 exponent_matrix，表示 x^1 y^2 z^3
    _exponent_matrix = Matrix([[1], [2], [3]])

    # 计算微分
    result = differential.calculate(_exponent_matrix)

    # 打印结果
    print("微分计算结果: ", result)
