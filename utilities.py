from sympy import Matrix as _Matrix, pprint
from sortedcontainers import SortedList
from math import ceil, sqrt


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


class Scalar:
    """
    In our setting, scalars are elements in a fintie prime Field of order (and characteristic) c
    """
    def __init__(self, c: int, val):
        self.c = c
        if val is not None:
            self.val = val % self.c
        else:
            self.val = None

    # Basic Functionalities

    def __eq__(self, other):
        return self.c == other.c and (self.val == other.val)

    def __str__(self):
        return str(self.val) + " (mod " + str(self.c) + ")"

    # Modifiers
    def increase_by(self, s: 'Scalar'):
        self.val = (self.val + s.val) % self.c

    def multiply_by(self, s: 'Scalar'):
        self.val = (self.val * s.val) % self.c

    def take_inverse(self):
        self.val = pow(self.val, -1, self.c)

    def update(self, val):
        self.val = val % self.c

    # Operations that yields new objects
    def __add__(self, other: 'Scalar'):
        return Scalar(self.c, self.val + other.val)

    def __sub__(self, other):
        return Scalar(self.c, self.val - other.val)

    def __mul__(self, other):
        return Scalar(self.c, self.val * other.val)

    def get_inverse(self):
        return Scalar(self.c, pow(self.val, -1, self.c))

    def __truediv__(self, other: 'Scalar'):
        return Scalar(self.c, self.val * pow(other.val, -1, self.c))

    # Static method to produce empty scalar wrapper
    @staticmethod
    def get_empty_scalar(c: int):
        return Scalar(c, None)


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
    def double_reduction(m1: "Matrix", m2: "Matrix" = None):
        if m2 is None:
            m2 = Matrix([[]]*m1.rows)
        output1 = []
        output2 = []
        combined = m1.copy().row_join(m2)
        for i in combined.rref()[1]:
            if i < m1.cols:
                output1.append(i)
            else:
                output2.append(i - m1.cols)
        return output1, output2

    def col_spans(self, vec: 'Vector') -> bool:
        if self.rows != vec.rows:
            raise ValueError("Number of rows do not match.")

        augmented_matrix = self.row_join(vec)

        # Compare ranks
        return self.rank() == augmented_matrix.rank()


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
        return scaled_res[0] // d, scaled_res[1] //d

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


if __name__ == "__main__":
    tb = Vector([
        [1, 3, 8, 5],
        [2, 4, 6, 2]
    ])
    tv = Vector([[100], [100]])
    for vec in convex_integral_combinations(tb, tv):
        assert tb * vec == tv
        for c in vec:
            assert c >= 0
        pprint(vec)
