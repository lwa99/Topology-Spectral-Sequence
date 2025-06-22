from sympy.abc import symbols
from sympy import Poly, QQ
from element import Bidegree

a = Bidegree([1, 2])
print(a)
b = Bidegree(a)
print(b)
