from scalar import Scalar
from abstractPage import AbstractPage
from element import Element


class Page(AbstractPage):
    def get_scalar(self, val) -> Scalar:
        return Scalar(self.c, val)


if __name__ == "__main__":
    p = Page(["x", "y"], [(1, 2), (3, 4)], 7)
    c_1 = p.get_scalar(14)
    print(c_1)
    e_1 = Element(p, [1, 2], p.get_scalar(3))
    print(e_1)
    e_2 = Element(p, [1, 2], p.get_scalar(-1))
    print((e_1 + e_2))
