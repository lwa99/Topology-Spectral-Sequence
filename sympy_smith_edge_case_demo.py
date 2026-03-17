"""
Minimal demonstration of a SymPy Smith normal form edge case on a field domain.

Run:
    python sympy_smith_edge_case_demo.py
"""

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp
from sympy.polys.domains import GF


def run_case(rows, domain):
    A = Matrix(rows)
    print(f"\nInput matrix A = {A.tolist()}, domain = {domain}")
    try:
        D, U, V = smith_normal_decomp(A, domain=domain)
        print("smith_normal_decomp succeeded.")
        print(f"D = {D.tolist()}")
        print(f"U = {U.tolist()}")
        print(f"V = {V.tolist()}")
    except Exception as exc:
        print("smith_normal_decomp raised an exception.")
        print(f"{type(exc).__name__}: {exc}")


if __name__ == "__main__":
    domain = GF(2)

    # Control case: same shape, no failure.
    run_case([[1, 0, 0]], domain)

    # Edge case observed in this project:
    # Expected to be easy over a field, but raises NotInvertible("zero divisor") in SymPy.
    run_case([[0, 0, 1]], domain)
