from __future__ import annotations
from sympy import Matrix as SMatrix, Poly as _Poly
from itertools import product
from warnings import warn


class Prime:
    """
    A static class used to handle prime number computation efficiently.
    """
    prime_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
                  89, 97, 101, 103, 107, 109, 113]

    @classmethod
    def is_prime(cls, n):
        from math import ceil, sqrt
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
                    cls.prime_list.append(d)
                d += 2
            return cls.prime_list[:]


class Matrix(SMatrix):
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
        non_empty_args = []
        for m in args:
            if m.rows != 0:
                non_empty_args.append(m)
        if len(non_empty_args) == 0:
            res: list[list | Matrix] = [[]]*len(args)
            return tuple(res + [Matrix()])

        rref, pivots = Matrix.hstack(*non_empty_args, Matrix.eye(non_empty_args[0].rows)).rref()
        acc_col = [0]
        res = []

        i = 0
        flag = False
        for j in pivots:
            while j >= acc_col[-1]:
                if i >= len(args):
                    flag = True
                    break
                acc_col.append(acc_col[-1] + args[i].cols)
                i += 1
                res.append([])
            if flag:
                break
            res[-1].append(j - acc_col[-2])
        res.extend([[]] * (len(args) - i))

        total_cols = sum([m.cols for m in non_empty_args])
        return tuple(res + [rref[:, total_cols:]])

    def __getitem__(self, item):
        """
        In case that we use a single index, make sure that the matrix is a row or column vector
        """
        try:
            iter(item)
            return super().__getitem__(item)
        except TypeError:
            if self.cols != 1 and self.rows != 1:
                warn("single indexing applied on non-vector matrix.")
            return super().__getitem__(item)

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


class Poly(_Poly):
    def __str__(self):
        return str(self.as_expr())

    __repr__ = __str__


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


def convex_integral_combinations(b: Matrix, v: Vector) -> tuple[tuple, ...]:
    """
    This function is devoted to solve the following problem:

    Let $b$ be a length n (>0) collection of 2-dimensional vectors and let $v$ be a specific 2-dimensional vector.
    Find all combinations of vectors in $b$ with non-negative integer coefficients that can sum up to $v$.
    We assume that the second component $v$ and those of vectors in $b$ are non-negative. Also, for vectors in b,
    we assume that when the second component is zero, their first component must be 0.

    Solution:
    Step 1: If $b$ contains only 1 vector $u$, test if $v$ is a multiple of $u$.
    Step 2: If $b$ contains exactly 2 or more vectors but all of them are dependent. Find all combinations using
    brute force
    Step 3: If there are at least 2 pivots.
        Step 3.1 If there are exactly 2 vectors: solve the linear system and check if the coefficients work out.
        Step 3.2 Calculate the bounds of coefficients corresponding to the free columns.
        (The bounds exists because the second components are non-negative).
        Step 3.3 Traverse through all linear combinations of the free columns within the bounds and check each case if
            the coefficients work out.
    """
    n = b.cols
    assert n > 0
    assert v[1] >= 0
    for i in range(n):
        assert b[1, i] >= 0

    # Step 1: Only one column.
    if n == 1:
        if b[1, 0] > 0:
            factor = v[1] // b[1, 0]
            if b.col(0) * factor == v:
                return ((factor,),)
            return ()
        else:  # b[1, 0] = 0
            factor = v[0] // b[0, 0]
            if b.col(0) * factor == v and factor >= 0:
                return ((factor,),)
            return ()

    res: list[tuple, ...] = []
    pivots = b.rref()[1]

    # Step 2: All columns dependent.
    if pivots == (0, ):
        if b[1, 0] == 0:  # Then all columns have the form (x, 0), x > 0
            if v[1] != 0 or v[0] < 0:
                return ()
            elif v[0] == 0:
                return (tuple([0]*n),)
            bounds = [v[0] // b[0, i] for i in range(n)]
        else:
            bounds = [v[1] // b[1, i] for i in range(n)]
        free_part_config = [0] * n

        try:
            while True:
                if b * Vector(free_part_config) == v:
                    res.append(tuple(free_part_config))
                free_part_config = _next_config(free_part_config, bounds)
        except StopIteration:
            return tuple(res)

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
            return ()
        return (config,)

    # Step 3.2 Calculate the bounds.
    bounds = [-1] * n
    skipped_index = []
    for j in range(n):
        if b[1, j] > 0:
            bounds[j] = v[1] // b[1, j]
            if b[0, j] > 0 and v[0] // b[0, j] < bounds[j]:
                bounds[j] = v[0] // b[0, j]
            continue
        skipped_index.append(j)

    # Now, we are left with those columns with 0 in the second component
    if len(skipped_index) > 0:
        cap = v[0]
        for j in range(n):
            if j in skipped_index:
                continue
            # collect negative first components of columns that were not skipped
            if b[0, j] < 0:
                cap -= b[0, j] * bounds[j]

        for j in skipped_index:
            if b[0, j] <= 0:
                raise ValueError  # The first component must be positive when the second is 0
            bounds[j] = cap // b[0, j]
            if bounds[j] < 0:
                return ()  # the cap and the first grade have different sign. There is no valid output

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
                res.append(tuple(config))
            free_part_config = _next_config(free_part_config, bounds)
    except StopIteration:
        return tuple(res)


def condition_fn(coefs, degree):
    # 后面需要什么条件自己再替换
    return True


def degree_generator(basis_matrix, c):
    basis_matrix = list(zip(*basis_matrix))

    n = len(basis_matrix)
    k = len(basis_matrix[0]) if n > 0 else 0

    for coefs in product(range(c), repeat=n):
        degree = [0] * k
        for i in range(n):
            for j in range(k):
                degree[j] += coefs[i] * basis_matrix[i][j]

        if not condition_fn(coefs, degree):
            continue

        yield Vector(degree)


if __name__ == "__main__":
    a = Matrix([[1, 2, 3], [4, 5, 6]])
    print(a[2])
