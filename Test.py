from utilities import Matrix, Vector, Polynomial, monomial
from sympy import pprint
from spectral_sequence import SpectralSequence
from page import Page
from element import HomoElem, Bigrade

ss = SpectralSequence(
    ["x", "y", "z"],
    Matrix([
        [7, 3, 0],
        [1, 0, 2]
    ]),
    97
)

x, y, z = ss.def_generators()

ss.kill(x**2, y**4, z**2)
# ss.add_relation(ss([2, 0, 0]))
# ss.add_relation(ss([0, 4, 0]))
# ss.add_relation(ss([0, 0, 2]))

Z = monomial([0, 0, 0])

print(ss.get_abs_basis(Bigrade([7, 3])))
p_1 = Page(
    ss, 1,
    {
        x: 0*Z,
        y: 0*Z,
        z: 0*Z
    },
    (1, -1)
)

p_2 = Page(
    ss, 2,
    {
        x: y**3
    },
    (2, -1)
)

p_3 = Page(
    ss, 3,
    {
        z: y
    },
    (3, -2)
)

m = p_3.get_module(Bigrade([9, 2]))

print(m.dim)
print(m.ker_basis)
print(ss.get_abs_basis(Bigrade([9, 2])))
for j in range(m.sp_basis.cols):
    print(m.sp_basis.col(j))
    print(HomoElem(m.page, abs_coordinate=m.sp_basis.col(j), abs_bigrade=m.bigrade))

