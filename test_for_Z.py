from typing import List, Optional, Tuple
from sympy import Matrix, zeros, ZZ
from sympy.matrices.normalforms import smith_normal_decomp
import sympy

print(sympy.__version__)

def check_well_defined(
    m: int,
    relations: List[List[int]],           # r_i in Z^m: generators of the submodule we mod out
    domain_vectors: List[List[int]],      # v_j in Z^m: points where you prescribed images
    images: List[List[int]]               # y_j in Z^n: images in Z^n (same length as domain_vectors)
) -> Tuple[bool, Optional[Matrix]]:
    """
    Decide if a homomorphism \bar{phi}: Z^m / <relations> -> Z^n sending [v_j] -> y_j is well-defined.
    If yes, also produce one lifting Phi: Z^m -> Z^n as an integer matrix A (n x m) with:
        A * r_i = 0    and    A * v_j = y_j.

    Inputs:
        m: int, dimension of the free domain Z^m.
        relations: list of r_i (each length m), generators of the submodule we quotient by.
        domain_vectors: list of v_j (each length m), the elements whose images are specified.
        images: list of y_j (each length n), images of v_j in Z^n.

    Returns:
        (is_well_defined, A) where:
            is_well_defined: bool
            A: sympy Matrix (n x m) with integer entries if True, else None
    """
    # basic shape checks
    if len(domain_vectors) != len(images):
        raise ValueError("domain_vectors and images must have the same length.")
    if any(len(r) != m for r in relations):
        raise ValueError("Each relation r_i must have length m.")
    if any(len(v) != m for v in domain_vectors):
        raise ValueError("Each domain vector v_j must have length m.")
    if len(images) == 0:
        # No constraints on images => a map exists (e.g., the zero map). Return zero witness.
        return True, zeros(0, m)  # n is unknown; with no images we can't infer n; returning empty witness.

    n = len(images[0])
    if any(len(y) != n for y in images):
        raise ValueError("All images y_j must have the same length n.")

    # Build U = [R | V] as an m x t integer matrix (columns are constraints)
    cols = []
    for r in relations:
        cols.append(Matrix(m, 1, r))
    for v in domain_vectors:
        cols.append(Matrix(m, 1, v))
    U = Matrix.hstack(*cols) if cols else Matrix.zeros(m, 0)  # m x t, where t = p+s

    # Build B = [0 | Y] as an n x t integer matrix matching the columns of U
    p = len(relations)
    s = len(domain_vectors)
    if U.shape[1] != p + s:
        raise RuntimeError("Internal dimension mismatch.")
    if p + s == 0:
        # Truly no constraints; any map works. Return zero matrix as witness.
        return True, Matrix.zeros(n, m)

    B_left = Matrix.zeros(n, p)
    B_right = Matrix.hstack(*[Matrix(n, 1, y) for y in images]) if s > 0 else Matrix.zeros(n, 0)
    B = Matrix.hstack(B_left, B_right)  # n x (p+s)

    # We need to solve for A (n x m) with A * U = B over Z.
    # This decouples by rows: for each k, find a row vector x_k (1 x m) s.t. x_k * U = b_k (1 x t).
    # Transpose to standard Ax=c form over Z:  (U^T) * x_k^T = b_k^T
    A_constraint = U.T  # shape (t x m)

    # SNF: P * A_constraint * Q = D  (over Z), with P,Q unimodular, D diagonal (rectangular)
    _snf = smith_normal_decomp(A_constraint, domain=ZZ)
    print(_snf)
    D, P, Q = _snf[:3]

    diag_len = min(D.rows, D.cols)
    diag_vals = [int(D[i, i]) for i in range(diag_len) if D[i, i] != 0]
    rank = len(diag_vals)

    witness_rows = []

    for row_idx in range(n):
        b_row = B[row_idx, :]                 # 1 x t
        c = Matrix(b_row).T                   # t x 1
        Pc = P * c                            # transform RHS to SNF coordinates

        # Solvability conditions:
        # For i < rank: (Pc)_i must be divisible by d_i
        # For i >= rank: (Pc)_i must be zero
        ok = True
        for i in range(rank):
            di = int(D[i, i])
            if di == 0:
                continue
            if int(Pc[i, 0]) % di != 0:
                ok = False
                break
        if ok:
            for i in range(rank, D.rows):
                if int(Pc[i, 0]) != 0:
                    ok = False
                    break
        if not ok:
            return False, None

        # Construct one integer solution:
        # D z = Pc, with z length m; set free variables to 0
        z = Matrix.zeros(Q.cols, 1)  # m x 1
        for i in range(rank):
            di = int(D[i, i])
            z[i, 0] = int(Pc[i, 0]) // di
        # x = Q * z solves A_constraint * x = c
        x = Q * z    # m x 1
        witness_rows.append(x.T)  # 1 x m

    # Stack rows to form A (n x m)
    A_witness = Matrix.vstack(*witness_rows) if witness_rows else Matrix.zeros(n, m)
    # Optional sanity check:
    assert (A_witness * U - B) == Matrix.zeros(n, U.shape[1])

    return True, A_witness


# --------- Tiny usage demo ---------
if __name__ == "__main__":
    # Example: m=2, quotient by <(1,-1)>, prescribe images of e1 and e2 consistently
    m = 2
    relations = [[1, -1]]               # we mod out by the submodule generated by (1, -1)
    domain_vectors = [[1, 0], [0, 1]]   # e1, e2
    images = [[1, 1], [1, 1]]           # want [e1] and [e2] to go to the same in Z^2
    ok, A = check_well_defined(m, relations, domain_vectors, images)
    print("Well-defined?", ok)
    if ok:
        print("One lifting A (rows act on Z^m):")
        print(A)
        # Check: A*(1,-1)^T == 0; and A*e1 == A*e2 == (1,1)^T

    # Inconsistent example: require e1 -> (1,1), e2 -> (0,0) while (1,-1) must map to 0
    images_bad = [[1, 1], [0, 0]]
    ok2, A2 = check_well_defined(m, relations, domain_vectors, images_bad)
    print("Well-defined (inconsistent)?", ok2)


