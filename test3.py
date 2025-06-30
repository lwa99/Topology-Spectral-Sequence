from sympy.abc import symbols
from sympy import GF, QQ, Poly

a, t = symbols(["a", "t"])
p = Poly(t, [t, a], domain=QQ)
print(p)

