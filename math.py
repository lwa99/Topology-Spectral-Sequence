"""
Implements Finite Prime Field
"""

from warnings import warn


class PField:
    debug = False

    @staticmethod
    def _isprime(n) -> bool:
        if n < 2:
            return False
        i = 2
        while i * i <= n:
            if n % i == 0:
                return False
            i += 1
        return True

    def _check_range(self, *args):
        if PField.debug:
            for i in args:
                if i >= self.n:
                    warn("PField: Input Element Out of Range")

    def __init__(self, n: int):
        if PField.debug and not PField._isprime(n):
            raise RuntimeError("PField: Not Prime Order")
        self.n = n

    def add(self, e_1: int, e_2: int):
        self._check_range(e_1, e_2)
        return (e_1 + e_2) % self.n

    def add_inv(self, e: int):
        self._check_range(e)
        return (-e) % self.n

    def mul(self, e_1: int, e_2: int):
        self._check_range(e_1, e_2)
        return (e_1 * e_2) % self.n

    def mul_inv(self, e: int):
        self._check_range(e)
        return pow(e, -1, self.n)


class BiGradedAlgebra:
    class HTerm:
        """
        Homogeneous Term

        Every homogeneous term is a function that maps
        """

    class HSpace:
        """
        Space of Homogeneous Terms
        """

    def __init__(self, basis_name: list[str], basis_bigrades: list[tuple[int, int]]):
        self.basis_name = basis_name
        self.basis_bigrades = basis_bigrades
        self.h_spaces: list[list[BiGradedAlgebra.HSpace]] = []
        # This 2-D array is ordered in both direction to facilitate further search.
        # TODO: It might be a better practice to define a class for this

    def calculate_bigrade(self):
        pass


if __name__ == "__main__":
    F = PField(7)
    print(F.add_inv(4))
    print(F.mul_inv(4))
