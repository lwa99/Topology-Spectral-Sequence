from sympy.abc import symbols
from sympy import GF, QQ
from spectral_sequence import SpectralSequence
from element import HomoElem, Bigrade

x, y, z = symbols(["x", "y", "z"])
ss = SpectralSequence(
    QQ,
    [x, y, z],
    [[7, 3, 0],
     [1, 0, 2]],
    [[1, 0],
     [-1, 1]]
)
# Note: [[1, 0],
#       [-1, 1]]
# means that the bigrade of d_n is
#       [[1,  0], *  [n  = [n
#        [-1, 1]]     1] =  -n + 1]

ss.kill(x**2, y**4, z**2)
p_1 = ss.add_page({x: 0, y: 0, z: 0})
p_2 = ss.add_page({x: y**3})
p_3 = ss.add_page({z: y})
p_4 = ss.add_page()
m = p_4[10, 1]

print(m.dim)
print(m.ker_basis)
for j in range(m.sp_basis.cols):
    print(HomoElem(m.page, abs_coordinate=m.sp_basis.col(j), abs_bigrade=m.bigrade))
