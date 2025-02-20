from scalar import Scalar
from abstractPage import AbstractPage
from element import Element
from numpy import ndarray, array


class Page(AbstractPage):
    def get_scalar(self, val) -> Scalar:
        return Scalar(self.c, val)


if __name__ == "__main__":
    try:
        assert False
    except AssertionError:
        print("Debug Mode: ON")

    p = Page(["x", "y"], array([[1, 2], [3, 4]]), 7)
    print(p.get_scalar(9))
    e_1 = Element(p, array([1, 2]), p.get_scalar(3))
    print(e_1)
    e_2 = Element(p, array([1, 2]), p.get_scalar(5))
    print((e_1 + e_2))
