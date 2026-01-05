# from sympy.abc import symbols
# from sympy import GF, QQ, Poly
#
# a, t = symbols(["a", "t"])
# p = Poly(t, [t, a], domain=QQ)
# print(p)

from sympy.polys.matrices.domainscalar import DomainScalar
from sympy import ZZ
a = DomainScalar(-1, ZZ)
print(a)
