from warnings import warn
class PField:
    """
    Fintie prime Field of characteristic (and order) n
    """
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