from tools import Matrix, pprint
from abstractPage import AbstractPage


def convex_integral_combinations():
    """
    This function is devoted to solve the following problem:

    Let $b$ be a length n (>0) collection of 2-dimensional vectors and let $v$ be a specific 2-dimensional vector.
    Find all combinations of vectors in $b$ with positive integer coefficients that can sum up to $v$.
    We assume that the first components of vectors in $b$ and both components of $v$ are non-negative.

    Solution:
    Step 1: If $b$ contains only 1 vector $u$, test if $v$ is a multiple of $u$.
    Step 2: If $b$ contains exactly 2 vectors $u_1, u_2$.
        Step 2.1 If $u_1, u_2$ are linearly dependent (decide by calculating determinant), use step 1
        Step 2.2 If they are linearly independent, solve the linear system and check if the coefficients work out.
    Step 3: If $b$ contains more than 2 vectors, $u_1, dots, u_n$
        Step 3.1 Scan through the list to see if there are two independent columns, if no, use step 1
        Step 3.2 Call the two independent columns "free". Calculate the bounds of coefficients corresponding to the
            free columns. (The bounds exists because the first components are positive).
        Step 3.3 Traverse through all linear combinations of the free columns within the bounds and check each case if
            the coefficients work out.
    """
    pass


def reduce_basis_information():
    """

    """
    pass


class Module:
    def __init__(self, page: AbstractPage, bigrade: Matrix):
        # initialization should only be called by page.getModule
        self.page = page
        self.page.calculated_modules.append(self)
        assert bigrade.shape == (2, 1)
        self.bigrade = bigrade
        self.basis = []
        # Start Generating Basis

        # 1 Only one generator
        if self.page.gen_num == 1:
            exponent = self.bigrade[0] // self.page.generator_bigrades[0, 0]
            if self.page.generator_bigrades[:0] * exponent == self.bigrade:
                self.basis.append(Matrix([exponent]))
                return

        # 2: Two or more generators
        # In current version, we assume that the first two bigrades are linearly independent.
        adj_A = Matrix([
            [self.page.generator_bigrades[1, 1], -self.page.generator_bigrades[0, 1]],
            [-self.page.generator_bigrades[1, 0], self.page.generator_bigrades[0, 0]]
        ])
        d = adj_A.det()  # det(adj_A) = det(A) in this case
        assert d != 0

        def check_to_add(_b: Matrix, tail: Matrix = None):
            # We want to verify that elements of A^{-1}b are positive integers. The method is to verify that
            # element in adj(A)b are positive multiples of der(A)
            scaled_first_two = adj_A * _b
            for i in scaled_first_two:
                if i * d < 0 or i % d != 0:
                    break
            else:
                if tail is None:
                    self.basis.append(scaled_first_two / d)
                else:
                    self.basis.append((scaled_first_two / d).col_join(tail))

        # 2.1 Exactly two generators
        if self.page.gen_num == 2:
            check_to_add(self.bigrade)
            return

        # 2.2 More than two generators
        bounds = self._calculate_bounds()[:-2]  # we only care about the bounds of last n - 2 variables
        current_exp = [0] * (self.page.gen_num - 2)
        while True:
            try:
                b = self.bigrade - self.page.generator_bigrades[:, 2:] * Matrix(current_exp)
                check_to_add(b, Matrix(current_exp))
                current_exp = Module._next_exponent(current_exp, bounds)
            except StopIteration:
                break

    def _calculate_bounds(self):
        """
        This method finds the bounds of the exponents of each variable
        """

        bg = self.page.generator_bigrades
        bounds = [-1] * self.page.gen_num
        skipped_index = []
        for j in range(self.page.gen_num):
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
        return bounds

    @staticmethod
    def _next_exponent(cur_exponent: list, bounds: list):
        # Find the last index that is not at the bound
        critical_index = len(bounds) - 1
        while cur_exponent[critical_index] == bounds[critical_index]:
            critical_index -= 1
            if critical_index == -1:
                raise StopIteration
        output = cur_exponent[:critical_index + 1]
        output[critical_index] += 1
        output.extend([0] * (len(bounds) - critical_index - 1))
        return output

    def get_basis(self):
        return self.basis


if __name__ == "__main__":
    # 生成测试 Page
    generators = ["x", "y", "z", "w"]
    generator_bigrades = Matrix([[1, 0, 4, 2], [4, -3, 6, 7]])  # 2x2 矩阵，每列是一个 generator 的 bigrade
    c = 7  # 质数特征
    print("Generator bigrades:")
    pprint(generator_bigrades)

    _page = AbstractPage(generators, generator_bigrades, c, 1)

    # 设定要测试的 bigrade
    _bigrade = Matrix([1000, 1000])
    print("Target:")
    pprint(_bigrade)

    # 初始化并计算 basis
    initializer = Module(_page, _bigrade)
    _basis = initializer.get_basis()

    # 打印 basis 结果
    print("Computed Basis:")
    for monomial in _basis:
        print(monomial)
