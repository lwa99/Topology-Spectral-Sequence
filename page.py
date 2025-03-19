from scalar import Scalar
from abstractPage import AbstractPage
from element import Element
from module import Module
from tools import Matrix, pprint, enumerate_exponents


class Page(AbstractPage):
    def set_differential(self):
        """
        Set the information of the differential
        """

        raise NotImplementedError

    def get_scalar(self, val) -> Scalar:
        return Scalar(self.c, val)

    def get_module(self, bigrade: Matrix):
        for module in self.calculated_modules.keys():
            if bigrade == module:
                return module
            return Module(self, bigrade)

    def _calculate_next_page_module(self):
        """
        Calculate the module of the next page at a give bigrade
        """
        # Step 1: Locate the "previous" module and the "next" module



        # Step 2: Calculate the matrix representations of differential


        # Step 3: Return kernel / image
        raise NotImplementedError

    def _compare_elements(self, e1, e2):
        raise NotImplementedError


    def compute_differential_matrix(self, domain_bigrade: Matrix) -> np.ndarray:
        """
        对于给定的 domain bigrade，利用当前设置好的 differential，
        计算 differential 在该向量空间上的矩阵表达。
    
        流程：
        1. 枚举出满足 self.generator_bigrades * exponent = domain_bigrade 的所有指数（即 domain basis）。
        2. 对于每个基元（以 Element 表示，系数取 1），计算其 differential。
        3. 根据第一个非零 differential 的 bigrade 确定 codomain basis。
        4. 将 differential 结果展开成 codomain basis 的线性组合，构造出矩阵表示。
        """
        # 1. 枚举 domain basis
        domain_exponents = enumerate_exponents(self, domain_bigrade)
        if not domain_exponents:
            print("未找到符合 domain bigrade 的基。")
            return np.array([])
        domain_basis = [Element(self, exp, self.get_scalar(1)) for exp in domain_exponents]

        # 2. 对每个 domain basis 元素应用 differential
        differential_images = [self.differential.calculate(exp) for exp in domain_exponents]

        # 3. 确定 codomain bigrade（取第一个非零 differential image 的 bigrade）
        target_bigrade = None
        for img in differential_images:
            if len(img.coefMap) > 0 and img.bigrade is not None and len(img.bigrade) > 0:
                target_bigrade = img.bigrade
                break
        if target_bigrade is None:
            print("所有 differential 结果均为零。")
            return np.zeros((0, len(domain_exponents)), dtype=int)
    
        # 4. 枚举 codomain basis
        codomain_exponents = enumerate_exponents(self, target_bigrade)
        if not codomain_exponents:
            print("未找到符合 codomain bigrade 的基。")
            return np.array([])
    
        # 构造矩阵：行对应 codomain basis, 列对应 domain basis
        M = np.zeros((len(codomain_exponents), len(domain_exponents)), dtype=int)
        for j, img in enumerate(differential_images):
            for key, scalar in img.coefMap.items():
                found = False
                for i, exp in enumerate(codomain_exponents):
                    if key == exp:
                        M[i, j] = scalar.val
                        found = True
                        break
                if not found:
                    print(f"警告：differential 中的项 {key} 不在 codomain basis 内。")
        return M



if __name__ == "__main__":
    try:
        assert False
    except AssertionError:
        print("Debug Mode: ON")

    p = Page(["x", "y", "z"], Matrix([[1, 2, 3], [3, 4, 5]]), 7, 1)
    print("Generator Bigrades: (Displayed as Columns)")
    pprint(p.generator_bigrades)
    e_1 = Element(p, Matrix([1, 2, 0]), p.get_scalar(3))
    print("e_1: ", e_1)
    e_2 = Element(p, Matrix([3, 1, 0]), p.get_scalar(5))
    e_3 = e_1 * e_2
    print("e_3: ", e_3)
    print("e_3 bigrade:")
    pprint(e_3.bigrade)
