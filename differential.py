from tools import Matrix
from page import Page
from element import Element
from scalar import Scalar


class Differential:
    def __init__(self, page: Page):
        """
        :param page: The associated Page instance.
        """
        self.page = page

    def calculate(self, diff_matrix, exponent_matrix: Matrix) -> Element:
        """
        Compute the differential of a monomial.

        :param diff_matrix: A list [[d(x)], [d(y)], [d(z)]] where each entry is an Element.
        :param exponent_matrix: Matrix representing the monomial exponents.
        :return: Element representing the calculated differential.
        """
        assert exponent_matrix.shape == (self.page.gen_num, 1)  # Ensure correct shape
        assert len(diff_matrix) == self.page.gen_num  # Ensure diff_matrix is valid

        differential_result = Element(self.page)

        for i in range(self.page.gen_num):
            exponent = exponent_matrix[i, 0]  # Get the exponent for this generator
            if exponent > 0:
                new_exponent = exponent_matrix.copy()
                new_exponent[i, 0] -= 1  # Reduce power by 1

                # Convert scalar coefficient to Element
                coefficient_scalar = self.page.get_scalar(exponent)  # Scalar 类型
                coefficient_element = Element(
                    self.page, 
                    Matrix.zeros(self.page.gen_num, 1),  # ✅ 确保是 (3,1) 矩阵
                    coefficient_scalar
                )  

                # Multiply the differential term by the scalar correctly
                differential_term = diff_matrix[i][0]  # Now accessing correctly
                term = differential_term * coefficient_element  # ✅ 直接用 Element 的乘法

                # 更新指数
                term.bigrade = self.page.generator_bigrades * new_exponent
                
                # ✅ 安全更新 coefMap，避免 KeyError
                if differential_term.bigrade in term.coefMap:
                    term.coefMap[new_exponent] = term.coefMap.pop(differential_term.bigrade)

                differential_result += term  # Sum up the differentials

        return differential_result




if __name__ == "__main__":
    # 创建 Page 实例
    generators = ["x", "y", "z"]
    generator_bigrades = Matrix([[1, 0, 3], [4, -3, 6]])  # 生成元的双次数
    c = 7  # 有限域特征数

    page = Page(generators, generator_bigrades, c)

    # 创建 Differential 实例
    differential = Differential(page)

    # 生成 d(x), d(y), d(z) 作为 Element
    dx = Element(page, Matrix([[1], [0], [0]]), Scalar(c, 1))  # d(x) 对应于 x^1 y^0 z^0
    dy = Element(page, Matrix([[0], [1], [0]]), Scalar(c, 1))  # d(y) 对应于 x^0 y^1 z^0
    dz = Element(page, Matrix([[0], [0], [1]]), Scalar(c, 1))  # d(z) 对应于 x^0 y^0 z^1

    # 这里使用普通的 Python 列表存储微分变量
    diff_matrix = [[dx], [dy], [dz]]

    # 选取 exponent_matrix，表示 x^1 y^2 z^3
    exponent_matrix = Matrix([[1], [2], [3]])

    # 计算微分
    result = differential.calculate(diff_matrix, exponent_matrix)

    # 打印结果
    print("微分计算结果: ", result)
