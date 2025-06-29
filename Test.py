from sympy.abc import symbols
from sympy import GF, QQ
from spectral_sequence import SpectralSequence
from element import HomoElem

a, t = symbols(["a", "t"])
ss = SpectralSequence(
    GF(3),
    [a, t],
    [[3, 0],  # Vertical
     [0, 2]],
    [[1, 0],
     [-1, 1]]
)

# Note: [[1, 0],
#       [-1, 1]]
# means that the bidegree of d_n is
#       [[1,  0], *  [n  = [n
#        [-1, 1]]     1] =  -n + 1]

ss.kill(a**2)
p_1 = ss.add_page({a: 0, t: 0})
p_2 = ss.add_page({a: 0, t: 0})
p_3 = ss.add_page({t: a})
p_4 = ss.add_page()
m = p_4[0, 6]

print(m.abs_dim)
print(m.dim)
print(m.ker_basis)
for j in range(m.sp_basis.cols):
    print(HomoElem(m.page, abs_coordinate=m.sp_basis.col(j), abs_bideg=m.bidegree).poly.as_expr())
