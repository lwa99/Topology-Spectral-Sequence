from __future__ import annotations
from sympy import Matrix as _Matrix
from sortedcontainers import SortedList, SortedDict
from math import ceil, sqrt
from copy import deepcopy
from itertools import product

class Prime:
    """
    A static class used to handle prime number computation efficiently.
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
    Matrix with dictionary order.
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
    def multi_reduction(*args: Matrix):
        rref, pivots = Matrix.hstack(*args, Matrix.eye(args[0].rows)).rref()
        acc_col = [0]
        res = []

        i = 0
        for j in pivots:
            if j >= acc_col[-1]:
                if i >= len(args):
                    break
                acc_col.append(acc_col[-1] + args[i].cols)
                i += 1
            res[i].append(j - acc_col[-2])

        total_cols = sum([m.cols for m in args])
        return tuple(res + [rref[:, total_cols:]])

    def col_spans(self, vec: 'Vector') -> bool:
        if self.rows != vec.rows:
            raise ValueError("Number of rows do not match.")

        augmented_matrix = self.row_join(vec)

        # Compare ranks
        return self.rank() == augmented_matrix.rank()

    def __repr__(self):
        return f"Matrix {self.tolist()}"

    def __str__(self):
        return self.__repr__()


class Vector(Matrix):
    def __new__(cls, *args, **kwargs):
        instance = super().__new__(cls, *args, **kwargs)
        if instance.cols != 1:
            return Matrix(*args, **kwargs)
        return instance

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__((key, 0))
        return super().__getitem__(key)

    def __repr__(self):
        return f"Vector {self.tolist()}"


class Exponent(Vector):
    def __repr__(self):
        return f"Exponent {self.tolist()}"


class Polynomial(SortedDict):
    base_field = None
    variables = None

    @classmethod
    def initiate(cls, base_field, variables):
        if Polynomial.base_field is not None:
            raise RuntimeWarning("Polynomial: Base Field Overridden")
        cls.base_field = base_field
        cls.variables = variables

    def __add__(self, other: Polynomial):
        assert self.base_field is not None, "Set the base field first."
        if len(self) == 0:
            return deepcopy(other)
        if len(other) == 0:
            return deepcopy(self)

        output = deepcopy(self)
        for exp, coef in other.items():
            if exp in output.keys():
                temp = output[exp] + coef
                if temp == self.base_field(0):
                    del output[exp]
                else:
                    output[exp] = temp
            else:
                output[exp] = coef
        return output

    def __mul__(self, other):
        assert self.base_field is not None, "Set the base field first."
        if isinstance(other, Polynomial):
            if len(self) == 0:
                return deepcopy(self)
            if len(other) == 0:
                return deepcopy(other)

            # From now on we may assume that both operands are non-zero
            output = Polynomial()
            for exp_1, coef_1 in self.items():
                for exp_2, coef_2 in other.items():
                    output[exp_1 + exp_2] = coef_1 * coef_2
        else:
            output = deepcopy(self)
            for exp, coef in output.items():
                output[exp] = coef * other
        return output

    def __rmul__(self, other):
        return self * other

    def __sub__(self, other):
        return self + self.base_field(-1) * other

    def __pow__(self, power, modulo=None):
        assert power >= 0
        if power == 0:
            return monomial([0]*len(self.variables))
        if power == 1:
            return self
        temp = (self ** (power // 2))
        if power % 2 == 0:
            return temp * temp
        else:
            return temp * temp * self

    def __repr__(self):
        output = ""
        if len(self) == 0:
            return "zero polynomial"
        for key, value in self.items():
            output += str(value.val)
            for i, exponent in enumerate(key):
                output += f"({self.variables[i]}^{exponent})"
            output += " + "
        output = output[:-2]
        return output

    def __str__(self):
        return self.__repr__()


def monomial(exp: list):
    return Polynomial({Exponent(exp): Polynomial.base_field(1)})


# _stored: dict[tuple[int, int], tuple[tuple[int]]] = {(0, 1): ((0,),)}
#
#
# def decompositions(_n: int, _k: int) -> tuple[tuple[int, ...]]:
#     if (_n, _k) in _stored.keys():
#         return _stored[(_n, _k)]
#     if _k == 1:
#         _stored[(_n, 1)] = ((_n,),)
#         return ((_n,),)
#     res_builder = []
#     for _i in range(0, _n+1):
#         for config in decompositions(_n - _i, _k-1):
#             res_builder.append(
#                 (_i, ) + config
#             )
#
#     _stored[(_n, _k)] = tuple(res_builder)
#     return tuple(res_builder)

def _next_config(cur_config: list, bounds: list):
    # Find the last index that is not at the bound
    critical_index = len(bounds) - 1
    while cur_config[critical_index] == bounds[critical_index]:
        critical_index -= 1
        if critical_index == -1:
            raise StopIteration
    output = cur_config[:critical_index + 1]
    output[critical_index] += 1
    output.extend([0] * (len(bounds) - critical_index - 1))
    return output


def convex_integral_combinations(b: Matrix, v: Vector) -> list[Vector]:
    """
    This function is devoted to solve the following problem:

    Let $b$ be a length n (>0) collection of 2-dimensional vectors and let $v$ be a specific 2-dimensional vector.
    Find all combinations of vectors in $b$ with non-negative integer coefficients that can sum up to $v$.
    We assume that the first components of vectors in $b$ and both components of $v$ are non-negative.

    Solution:
    Step 1: If $b$ contains only 1 vector $u$, test if $v$ is a multiple of $u$.
    Step 2: If $b$ contains exactly 2 or more vectors but all of them are dependent. Find all combinations using
    brute force
    Step 3: If there are at least 2 pivots.
        Step 3.1 If there are exactly 2 vectors: solve the linear system and check if the coefficients work out.
        Step 3.2 Calculate the bounds of coefficients corresponding to the free columns.
        (The bounds exists because the first components are positive).
        Step 3.3 Traverse through all linear combinations of the free columns within the bounds and check each case if
            the coefficients work out.
    """
    n = b.cols
    assert n > 0
    assert v[0] >= 0 and v[1] >= 0
    for i in range(n):
        assert b[0, i] >= 0

    # Step 1: Only one column.
    if n == 1:
        factor = v[0] // b[0, 0]
        if b.col(0) * factor == v:
            return [Vector([factor])]

    res: list[Vector, ...] = []
    pivots = b.rref()[1]

    # Step 2: All columns dependent.
    if pivots == (0, ):
        bounds = [v[0] // b[0, i] for i in range(n)]
        free_part_config = [0] * n

        try:
            while True:
                if b * Vector(free_part_config) == v:
                    res.append(Vector(free_part_config))
                free_part_config = _next_config(free_part_config, bounds)
        except StopIteration:
            return res

    # Step 3 At least two pivots
    p0 = pivots[0]
    p1 = pivots[1]

    adj_a = Matrix([
        [b[1, p1], -b[0, p1]],
        [-b[1, p0], b[0, p0]]
    ])
    d = adj_a.det()

    def check(target: Vector) -> tuple[int, int] | None:
        scaled_res = adj_a * target
        for _i in scaled_res:
            if _i * d < 0 or _i % d != 0:
                return None
        return scaled_res[0] // d, scaled_res[1] // d

    # Step 3.1 Exactly two columns
    if n == 2:
        config = check(v)
        if config is None:
            return []
        return [Vector(config)]

    # Step 3.2 Calculate the bounds.
    bounds = [-1] * n
    skipped_index = []
    for j in range(n):
        if b[0, j] > 0:
            bounds[j] = v[0] // b[0, j]
            if b[1, j] > 0 and v[1] // b[1, j] < bounds[j]:
                bounds[j] = v[1] // b[1, j]
            continue
        skipped_index.append(j)

    # Now, we are left with those columns with 0 in the first row
    if len(skipped_index) > 0:
        cap = v[1]
        for j in range(n):
            if j in skipped_index:
                continue
            # collect second grades that have different sign from the first skipped second grade
            if b[1, j] * b[1, skipped_index[0]] < 0:
                cap -= b[1, j] * bounds[j]

        for j in skipped_index:
            if b[1, j] * b[1, skipped_index[0]] <= 0:
                raise ValueError  # if different signs present, there are infinitely many possibilities
            bounds[j] = cap // b[1, j]
            if bounds[j] < 0:
                return []  # the cap and the second grade have different sign. There is no valid output

    # Step 3.3 Traverse all possible combinations.
    del bounds[p1]
    del bounds[p0]
    free_part_config = [0] * (n-2)
    free_part_idx = list(range(n))
    del free_part_idx[p1]
    del free_part_idx[p0]

    try:
        while True:
            fix_part_config = check(v - b[:, free_part_idx] * Vector(free_part_config))
            if fix_part_config is not None:
                config = free_part_config.copy()
                config.insert(p0, fix_part_config[0])
                config.insert(p1, fix_part_config[1])
                res.append(Vector(config))
            free_part_config = _next_config(free_part_config, bounds)
    except StopIteration:
        return res
    
def condition_fn(coeffs, degree):
    # 后面需要什么条件自己再替换
    return True

def degree_generator(basis_matrix, c):
    basis_matrix = list(zip(*basis_matrix))

    n = len(basis_matrix)
    k = len(basis_matrix[0]) if n > 0 else 0

    for coeffs in product(range(c), repeat=n):
        degree = [0] * k
        for i in range(n):
            for j in range(k):
                degree[j] += coeffs[i] * basis_matrix[i][j]

        if not condition_fn(coeffs, degree):
            continue

        yield Vector(degree)


if __name__ == "__main__":

    basis = [
        [1, 2],
        [-1, 3],
        [0, 0]
    ]
    c = 2

    gen = degree_generator(basis, c)

    print(next(gen))
    print(next(gen))
    print(next(gen))
    print(next(gen))

 



    # tb = Vector([
    #     [1, 3, 8, 5],
    #     [2, 4, 6, 2]
    # ])
    # tv = Vector([[100], [100]])
    # for t_vec in convex_integral_combinations(tb, tv):
    #     assert tb * t_vec == tv
    #     for t_c in t_vec:
    #         assert t_c >= 0
    #     pprint(t_vec)

    # from sympy import GF
    # F = GF(7)
    # Polynomial.set_base_field(F)
    # p_1 = monomial([1, 2])
    # p_2 = Polynomial({Exponent([1, 2]): F(6)})
    # print(p_1, p_2)
    # print(p_1 - p_2)

    tv = Vector([1, 2, 3])
    print(type(tv[:-1]))

