from sympy import Matrix as _Matrix, pprint
from sortedcontainers import SortedList
from math import ceil, sqrt


class Prime:
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
    __eq__ is rewritten to compare the whole array.
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
            res *= pow(p, self[i])
        return int(res)


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
