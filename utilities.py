from sympy import Matrix as _Matrix, pprint
from sortedcontainers import SortedList, SortedDict
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
    def double_reduction(m1: "Matrix", m2: "Matrix"):
        output1 = output2 = []
        combined = m1.copy().row_join(m2)
        for i in combined.rref()[1]:
            if i > m1.cols:
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
    def __init__(self, *args, **kwargs):
        super.__init__(*args, **kwargs)
        assert self.cols == 1


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


if __name__ == "__main__":
    a = Matrix([[1, 2], [2, 3], [3, 4]])
    b = Matrix([[1, 2], [2, 3], [3, 4]])
    print(a == b)
    c = Matrix([[1, 2], [2, 3], [3, 3]])
    print(c < a)
    print(c.col_spans(Matrix([[3], [5], [7]])))

    a = Matrix([[1, 2], [3, 4]])
    b = Matrix([4, 3])
    pprint(a)
    pprint(b)
