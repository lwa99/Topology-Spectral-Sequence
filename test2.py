from sympy.abc import symbols
from sympy import Poly, QQ
from element import Bigrade

a = Bigrade([1, 2])
print(a)
b = Bigrade(a)
print(b)
