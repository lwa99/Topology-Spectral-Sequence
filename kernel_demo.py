from sympy.polys.domains import ZZ

from AlgebraicCore import SNFMatrix
from utilities import Matrix


def show_kernel(A_data):
    A = SNFMatrix(A_data, domain=ZZ)
    K = A.get_kernel()
    print("A =", Matrix(A.tolist()))
    print("rank(A) =", A.r)
    print("Kernel basis (columns):\n", K)
    # Verify A * K == 0
    if K.cols == 0:
        print("A * K = [] (trivial kernel)")
    else:
        print("A * K =\n", A * K)
    print("-" * 40)


def main():
    # Example 1: 2x2 rank-1 matrix
    show_kernel([[1, 2], [2, 4]])

    # Example 2: 1x3 row vector, rank 1, kernel dimension 2
    show_kernel([[1, 2, 3]])

    # Example 3: full rank square matrix (trivial kernel)
    show_kernel([[2, 0], [0, 3]])


if __name__ == "__main__":
    main()

