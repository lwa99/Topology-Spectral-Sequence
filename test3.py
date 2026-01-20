# from sympy.abc import symbols
# from sympy import GF, QQ, Poly
#
# a, t = symbols(["a", "t"])
# p = Poly(t, [t, a], domain=QQ)
# print(p)

# from sympy.polys.matrices.domainscalar import DomainScalar
# from sympy import ZZ
# a = DomainScalar(-1, ZZ)
from matrices import DM
from sympy import ZZ
M = DM([[1, 0, 0], [0, 1, 0], [-1, -1, 1]], ZZ)
K = DM([[1, 0, 0], [0, 1, 0], [1, 1, 1]], ZZ)
print(M.inv_den())
print(M*K)
