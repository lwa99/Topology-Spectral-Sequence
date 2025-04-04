from utilities import Matrix, Vector, Polynomial, monomial
from sympy import pprint
from spectral_sequence import SpectralSequence
from page import Page
from element import HomoElem

ss = SpectralSequence(
    ["k", "y", "z"],
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

input_str = input("HAHA?")
for j, g in enumerate(ss.generators):
    exponent = [0] * len(ss.generators)
    exponent[j] = 1
    exec(f"{g} = monomial(exponent)")
s = monomial([0]*len(ss.generators))
temp_poly = None
exec("temp_poly =  " + input_str)
print(temp_poly)
