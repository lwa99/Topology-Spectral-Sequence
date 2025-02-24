from scalar import Scalar
from abstractPage import AbstractPage
from element import Element
from tools import Matrix, pprint


class Page(AbstractPage):
    def get_scalar(self, val) -> Scalar:
        return Scalar(self.c, val)


if __name__ == "__main__":
    try:
        assert False
    except AssertionError:
        print("Debug Mode: ON")

    p = Page(["x", "y", "z"], Matrix([[1, 2, 3], [3, 4, 5]]), 7)
    print("Generator Bigrades: (Displayed as Columns)")
    pprint(p.generator_bigrades)
    e_1 = Element(p, Matrix([1, 2, 0]), p.get_scalar(3))
    print("e_1: ", e_1)
    e_2 = Element(p, Matrix([3, 1, 0]), p.get_scalar(5))
    e_3 = e_1 * e_2
    print("e_3: ", e_3)
    print("e_3 bigrade:")
    pprint(e_3.bigrade)
