from sympy import GF, poly, Poly, symbols
from utilities import Vector

x, y = symbols("x, y")
my_dict = {
    (1, 2): 3,
    (4, 5): 4
}


# print(Poly(x**3 + 2*x*y - 1, x, y).terms())
my_poly = Poly.from_dict(my_dict, x, y)
print(my_poly)