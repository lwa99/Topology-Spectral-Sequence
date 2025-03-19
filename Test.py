from utilities import Matrix, Vector, Polynomial
from sympy import pprint
from spectral_sequence import SpectralSequence
from page import Page
from element import HomoPoly

ss = SpectralSequence(
    ["x", "y", "z"],
    Matrix([
        [7, 3, 0],
        [1, 0, 2]
    ]),
    97
)
ss.add_relation(ss([2, 0, 0]))
ss.add_relation(ss([0, 4, 0]))
ss.add_relation(ss([0, 0, 2]))

p = Page(ss, 1)


test_poly = ss([1, 3, 1]) * ss(2)

abs_bigrade, abs_coordinate = ss.get_absolute_info(test_poly)
abs_basis = ss.get_actual_basis(abs_bigrade)

print("______ Basis", abs_basis)
print("______ Absolute Info", abs_bigrade, abs_coordinate)
print("______ Kernel Basis", ss.get_ker_basis(abs_bigrade))

homo_poly = HomoPoly(p, test_poly)
print(homo_poly.coordinate)
