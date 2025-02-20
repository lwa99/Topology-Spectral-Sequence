import numpy as np


class AbstractPage:
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

    def __init__(self, generators: list[str], generator_bigrades: np.ndarray, c: int):
        """
        :param generators: Names of the generators
        :param generator_bigrades: Bigrades of the generators
        :param c: Order (and characteristic) of the base field: need to be prime
        """
        assert len(generators) == len(generator_bigrades)
        self.generators = generators
        self.generator_bigrades = generator_bigrades
        assert self._isprime(c)
        self.c = c

    @property
    def gen_num(self) -> int:
        return len(self.generators)

    def get_scalar(self, val):
        raise NotImplementedError
