from sympy.abc import symbols
from sympy import Poly, QQ


class MyPoly(Poly):
    def __str__(self):
        return "s"

    def __repr__(self):
        return "r"


x, y, z = symbols(["x", "y", "z"])


print(MyPoly(x+y, domain=QQ))
