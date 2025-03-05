from tools import Matrix, Prime


class AbstractPage:
    def __init__(self, generators: list[str], generator_bigrades: Matrix, c: int, page_num):
        """
        :param generators: Names of the generators
        :param generator_bigrades: Bigrades of the generators
        :param c: Order (and characteristic) of the base field: need to be prime
        """
        assert len(generators) == len(generator_bigrades) / 2
        self.generators = generators
        self.generator_bigrades: Matrix = generator_bigrades
        assert Prime.is_prime(c)
        self.c = c
        self.page_num = page_num
        self.calculated_modules = []

    @property
    def gen_num(self) -> int:
        return len(self.generators)

    def get_scalar(self, val):
        raise NotImplementedError
