from sympy import Matrix as _Matrix, pprint
from sortedcontainers import SortedList
from math import ceil, sqrt


class Prime:
    """
    Class that handle prime number computation efficiently.
    """
    prime_list = SortedList([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                             53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113])

    @classmethod
    def is_prime(cls, n):
        if n in cls.prime_list:
            return True
        elif n < cls.prime_list[-1]:
            return False
        else:
            bound = ceil(sqrt(n))
            for p in cls.prime_list:
                if n % p == 0:
                    return False
                if p > bound:
                    return True

            d = cls.prime_list[-1] + 2
            while d <= bound:
                if n % d == 0:
                    return False
                d += 2
            return True

    @classmethod
    def first_n_prime(cls, n):
        if n <= len(cls.prime_list):
            return cls.prime_list[0:n]
        else:
            d = cls.prime_list[-1] + 2
            while len(cls.prime_list) < n:
                if cls.is_prime(d):
                    cls.prime_list.add(d)
                d += 2
            return cls.prime_list[:]


class Matrix(_Matrix):
    """
    Ndarray with dictionary order.
    __eq__ is rewritten to compare the whole matrix.
    __hash__ calls the hash function for tuples.
    """

    def __lt__(self, other):
        # noinspection PyTypeChecker
        for i, n in enumerate(self):
            if n < other[i]:
                return True
        return False

    def __hash__(self):
        res = 1
        for i, p in enumerate(Prime.first_n_prime(len(self))):
            if self[i] >= 0:
                res *= pow(p, 2 * self[i])
            else:
                res *= pow(p, -2 * self[i] - 1)  # -1 to 1, -2 to 3, -3 to 5 ...
        return int(res)

    @staticmethod
    def double_reduction(m1: "Matrix", m2: "Matrix"):
        output1 = output2 = []
        combined = m1.copy().row_join(m2)
        for i in combined.rref()[1]:
            if i > m1.cols:
                output1.append(i)
            else:
                output2.append(i - m1.cols)
        return output1, output2


class Vector:
    """
    A vector must implement a method that returns its coordinate on the standard basis as a column vector.
    """
    @property
    def coordinate(self) -> Matrix:
        raise NotImplementedError


class VectorSpace:
    """
    A vector space stores a basis of itself and a kernel subspace.
    """
    def __init__(self, basis: Matrix, ker_basis: Matrix):
        self._basis = basis
        self._ker_basis = ker_basis

    def enlarge_to_basis(self, partial_basis: Matrix):
        """
        Give a set of vectors, we want to
        """
        pass

    def get_quotient_space(self, quotient_basis):
        pass

    def __contains__(self, item):
        raise NotImplementedError

    def get_coordinate(self, vec: Vector):
        """
        Return the coordinate of a vector on current basis.
        """
        if vec not in self:
            return None
        return self._basis.inv() * vec.coordinate


class LinearTransformation:
    """
    We store a linear transformation by its matrix representation. It can be defined by its image
    on a specific set of vectors.
    """

    def __init__(self, images):
        self.images = images

    def get_matrix(self, from_space, to_space):
        pass


def enumerate_exponents(page: Page, target_bigrade: Matrix):
    """
    枚举所有满足:
         page.generator_bigrades * exponent == target_bigrade
    的非负整数指数组合。返回值为一列表，每个元素为形如 Matrix([[e1], [e2], ... [en]]) 的指数列，
    其中 n = page.gen_num。
    """
    results = []
    n = page.gen_num
    # 设 generator_bigrades 的行数为 r, 则 target_bigrade 的形状应为 (r, 1)
    r = page.generator_bigrades.shape[0]

    def rec(i, current):
        if i == n:
            exp_matrix = Matrix([[x] for x in current])
            # 判断是否满足：generator_bigrades * exp_matrix == target_bigrade
            if page.generator_bigrades * exp_matrix == target_bigrade:
                results.append(exp_matrix)
            return
        # 对于每个生成元 i，计算一个粗略的上界：
        bounds = []
        for j in range(r):
            gen_val = page.generator_bigrades[j, i]
            if gen_val > 0:
                candidate = target_bigrade[j, 0] // gen_val
                bounds.append(candidate)
        # 如果在某些坐标上没有正贡献，则我们限定该生成元指数为 0
        max_exp = min(bounds) if bounds else 0
        for exp in range(max_exp + 1):
            rec(i + 1, current + [exp])
    rec(0, [])
    return results


if __name__ == "__main__":
    a = Matrix([[1, 2], [2, 3], [3, 4]])
    b = Matrix([[1, 2], [2, 3], [3, 4]])
    print(a == b)
    c = Matrix([[1, 2], [2, 3], [3, 3]])
    print(c < a)

    a = Matrix([[1, 2], [3, 4]])
    b = Matrix([4, 3])
    pprint(a)
    pprint(b)
