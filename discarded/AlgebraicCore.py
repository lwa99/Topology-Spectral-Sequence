from __future__ import annotations

from utilities import Matrix, convex_integral_combinations
from sympy.matrices.normalforms import smith_normal_decomp
from sympy.polys.domains.domain import Domain
from sympy.polys.polyerrors import ExactQuotientFailed
from sympy.polys.domains import ZZ, QQ
from sympy.matrices.matrices import MatrixBase
from typing import Optional, List, Dict, Tuple, TYPE_CHECKING
from sortedcontainers import SortedDict

if TYPE_CHECKING:
    from src.element import HomoElem, Bidegree
    from page import Page

_verify = True


def _is_field(domain: Domain) -> bool:
    """Return True when the domain advertises field semantics."""
    return bool(getattr(domain, "is_Field", False) or getattr(domain, "is_field", False))


def divides(a, d, domain: Domain):
    """Return True iff d divides a in the domain (fields always succeed when d != 0)."""
    d_dom = domain.convert(d)
    a_dom = domain.convert(a)
    if d_dom == domain.zero:
        # By convention 0 only divides 0.
        return a_dom == domain.zero

    # Fast path for fields: every nonzero element divides any element
    if _is_field(domain):
        return True

    try:  # try exact division
        domain.exquo(a_dom, d_dom)
        return True
    except ExactQuotientFailed:
        return False


def _hstack(domain, *args) -> "SNFMatrix":
    """
    Helper to horizontally stack column vectors/matrices and wrap as SNFMatrix.
    """
    return SNFMatrix(Matrix.hstack(*args), domain=domain)


def homo_elem_to_coordinate(homo_elem: HomoElem) -> Optional[Matrix]:
    """Convert HomoElem to its coordinate vector (or None if zero).

    Args:
        homo_elem: Element to convert

    Returns:
        Column vector coordinate, or None if element is zero
    """
    if homo_elem.isZero():
        return None
    return homo_elem.coordinate


def coordinate_to_homo_elem(page: Page, bidegree, coordinate: Matrix):
    """Convert coordinate vector and bidegree to HomoElem.

    Args:
        page: Page containing the element
        bidegree: Bidegree of the element
        coordinate: Column vector coordinate in absolute basis

    Returns:
        HomoElem with given bidegree and coordinate
    """
    if coordinate is None or coordinate.is_zero_matrix:
        from src.element import HomoElem
        return HomoElem(page, expr="0")

    from src.element import HomoElem
    return HomoElem(page, abs_bideg=bidegree, abs_coordinate=coordinate)


def _compute_quotient_basis(span_set: Matrix, zero_set: Matrix, domain) -> Matrix:
    """
    Compute basis for quotient span_set / zero_set.

    Returns the column vectors that form a basis for the quotient module.
    If zero_set spans a submodule of span_set, this returns an independent
    set that spans the quotient.

    Args:
        span_set: Columns are basis vectors of full module
        zero_set: Columns are basis vectors of submodule (relations)
        domain: ZZ or QQ

    Returns:
        Matrix whose columns form basis for quotient

    Example:
        >>> span = Matrix([[1, 0, 1], [0, 1, 1]])  # 3 generators in Z^2
        >>> zero = Matrix([[1], [1]])  # 1 relation
        >>> basis = _compute_quotient_basis(span, zero, QQ)
        >>> basis.shape == (2, 2)  # Dimension 3 - 1 = 2
    """
    from sympy import Matrix

    # Edge cases
    if span_set.cols == 0:
        return Matrix.zeros(span_set.rows, 0)

    if zero_set.cols == 0:
        # No relations - entire span is quotient
        return span_set

    if _is_field(domain):
        # Field case: use row reduction
        # Stack [zero_set | span_set] and reduce to find quotient basis
        combined = zero_set.row_join(span_set)
        rref, pivots = combined.rref()

        # Pivots in the span_set part (after zero_set columns) give quotient basis
        span_offset = zero_set.cols
        quotient_pivots = [p - span_offset for p in pivots if p >= span_offset]

        if quotient_pivots:
            return span_set[:, quotient_pivots]
        else:
            # All span columns are dependent on zero_set
            return Matrix.zeros(span_set.rows, 0)
    else:
        # PID case (ZZ): Use kernel-like approach
        # The quotient basis is found by:
        # 1. Find embedding: injective map zero_set → span_set
        # 2. Compute kernel of this map
        # For simplicity, use SNF on the augmented system

        # Build map: identify which columns of span_set contain the zero_set
        # For now, assume zero_set is a submatrix of span_set (columns)
        # Compute a complement basis

        # Use kernel approach: find basis for ker of the relation matrix
        # Build the "quotient map": embed zero_set into span_set coordinate system
        combined = zero_set.row_join(span_set)
        snf = SNFMatrix(combined, domain=domain)

        # The quotient basis corresponds to free variables in span_set
        # Use kernel of zero_set in span_set's coordinate space
        m = zero_set.cols
        n = span_set.cols

        # For each generator in span_set, check if it's independent mod zero_set
        quot_cols = []
        for j in range(n):
            col_j = span_set[:, j]
            # Check if col_j is in span of previously selected quotient columns
            if len(quot_cols) == 0:
                # First column always independent
                quot_cols.append(col_j)
            else:
                current_basis = Matrix.hstack(*quot_cols)
                snf_test = SNFMatrix(current_basis, domain=domain)
                if not snf_test.column_spans(col_j):
                    quot_cols.append(col_j)

        if quot_cols:
            return Matrix.hstack(*quot_cols)
        else:
            return Matrix.zeros(span_set.rows, 0)


class SNFMatrix(Matrix):
    """Matrix that caches SNF (or field rank/nullspace) for kernel/span queries."""

    def __new__(cls, *args, domain=ZZ, **kwargs) -> "SNFMatrix":
        """Build an SNFMatrix; use field rank/nullspace or PID Smith normal form."""
        if len(args) == 1 and isinstance(args[0], MatrixBase):
            # Build from its data instead of passing the matrix itself
            self = super().__new__(cls, args[0].tolist(), **kwargs)
        else:
            self = super().__new__(cls, *args, **kwargs)

        # Verify entries are coercible to domain
        if _verify:
            for entry in self:
                domain.convert(entry)

        self.domain = domain

        # Fast path when the domain is a field: avoid SNF and use rank/nullspace
        if _is_field(domain):
            # Over a field, rank is simply the matrix rank
            self.r = int(super(SNFMatrix, self).rank())  # ensure int
            return self

        # PID path: use Smith normal form
        if _verify:
            assert domain.is_PID
        self._D, self._U, self._V = smith_normal_decomp(self, domain=domain)

        # Compute rank from SNF diagonal: number of non-zero diagonal entries
        r = 0
        for i in range(min(self.rows, self.cols)):
            if self._D[i, i] != 0:
                r += 1
            else:
                break
        self.r = r
        return self

    @classmethod
    def from_columns(cls, domain: Domain, dim: int, columns) -> "SNFMatrix":
        """Build from column vectors; empty input yields a dim x 0 zero matrix."""
        columns = list(columns)
        if not columns:
            base = Matrix.zeros(dim, 0)
            return cls(base, domain=domain)

        for v in columns:
            if v.rows != dim or v.cols != 1:
                raise ValueError(
                    "All columns must be column vectors of shape (dim, 1); "
                    f"got shape ({v.rows}, {v.cols})."
                )
        base = Matrix.hstack(*columns)
        return cls(base, domain=domain)

    def get_kernel(self):
        """Return a matrix whose columns form a basis of ker(self) over the domain."""
        if _is_field(self.domain):
            # Over a field, compute kernel via nullspace
            basis = super(SNFMatrix, self).nullspace()
            if not basis:
                return Matrix.zeros(self.cols, 0)
            return Matrix.hstack(*basis)

        # Using U * A * V = D, the kernel of A is generated by the
        # columns of V corresponding to the zero diagonal entries of D.
        # If r is the rank, then indices [r, ..., n-1] span ker(A).
        if self.r == self.cols:
            # Trivial kernel: return an empty matrix with correct row size
            return Matrix.zeros(self.cols, 0)
        return self._V[:, self.r:]

    def get_columnspace(self):
        """Return a matrix whose columns form a basis of col(self) over the domain.

        For fields (QQ): Uses standard columnspace via Gaussian elimination.
        For PIDs (ZZ): Uses SNF - if U*A*V = D, then col(A) = U[:, :r] where r = rank.

        Returns: Matrix of shape (rows, rank) whose columns form a basis of column space.
        """
        if _is_field(self.domain):
            # Over a field, use standard columnspace
            basis = super(SNFMatrix, self).columnspace()
            if not basis:
                return Matrix.zeros(self.rows, 0)
            return Matrix.hstack(*basis)

        # Using U * A * V = D, the column space of A is generated by the
        # columns of U^{-1} * D corresponding to non-zero diagonal entries.
        # Preserve torsion by scaling with the diagonal of D.
        if self.r == 0:
            return Matrix.zeros(self.rows, 0)
        U_inv = self._U.inv()
        cols = []
        for i in range(self.r):
            col = Matrix(U_inv[:, i]) * self._D[i, i]
            if _verify:
                for j in range(col.rows):
                    col[j, 0] = self.domain.convert(col[j, 0])
            cols.append(col)
        return Matrix.hstack(*cols)

    def column_spans(self, b: Matrix) -> bool:
        """Return True iff b is in the column span (rank test for fields, SNF for PIDs)."""
        assert b.cols == 1 and b.rows == self.rows
        if _is_field(self.domain):
            # Over a field, b is in the column span iff rank(A|b) == rank(A)
            A = self
            return int(A.rank()) == int(Matrix.hstack(A, b).rank())

        bprime = self._U * b

        # 1) Zero conditions on the lower block: b'_{i} = 0 for i > r
        for i in range(self.r, self.rows):
            if self.domain.convert(bprime[i, 0]) != self.domain.zero:
                return False

        # 2) Divisibility conditions on the first r coordinates: d_i | b'_i
        for i in range(self.r):
            di = self._D[i, i]
            if not divides(bprime[i, 0], di, self.domain):
                return False
        return True

    @staticmethod
    def multi_reduction_pid(domain: Domain, *args: Matrix):
        """
        PID-aware analogue of Matrix.multi_reduction using SNF-based span checks.

        Returns pivot indices per input block (in original column order) and the
        inverse of the combined basis matrix. We avoid field RREF so torsion information
        is preserved.
        """
        non_empty_args = [m for m in args if m.rows != 0]
        if len(non_empty_args) == 0:
            res: list[list | Matrix] = [[]] * len(args)
            # Return empty pivots and identity matrix
            return tuple(res + [Matrix.eye(args[0].rows if args else 0)])

        rows = non_empty_args[0].rows
        pivots_per_block: list[list[int]] = [[] for _ in args]
        selected_cols: list[Matrix] = []
        basis_mat = SNFMatrix.from_columns(domain, rows, [])

        for block_idx, block in enumerate(args):
            if block.rows == 0:
                continue
            for j in range(block.cols):
                col = block.col(j)
                if basis_mat.cols > 0 and basis_mat.column_spans(col):
                    continue
                pivots_per_block[block_idx].append(j)
                selected_cols.append(col)
                basis_mat = SNFMatrix.from_columns(domain, rows, selected_cols)

        # Compute the inverse of the selected basis matrix
        if len(selected_cols) == 0:
            # No columns selected, return identity
            inv = Matrix.eye(rows)
        else:
            combined_basis = Matrix.hstack(*selected_cols)

            # The selected columns should form a linearly independent set
            # spanning the combined module. This should typically be square.
            if combined_basis.cols > combined_basis.rows:
                raise ValueError(
                    f"multi_reduction_pid produced more columns than rows: "
                    f"{combined_basis.rows}x{combined_basis.cols}. "
                    f"This indicates the input matrices have more independent columns than ambient dimension."
                )

            if combined_basis.cols == combined_basis.rows:
                # Square matrix: compute actual inverse using domain-appropriate method
                if _is_field(domain):
                    # Over fields: use standard Gaussian elimination-based inverse
                    inv = combined_basis.inv()
                else:
                    # Over PIDs: use SNF-based inverse for exactness
                    snf_basis = SNFMatrix(combined_basis, domain=domain)
                    if snf_basis.r != combined_basis.rows:
                        raise ValueError(
                            "Basis matrix is not invertible over PID: "
                            f"rank {snf_basis.r} != {combined_basis.rows}"
                        )
                    # Compute inverse using SNF: if U*A*V = D, then A^{-1} = V*D^{-1}*U
                    # For square matrix with full rank, D is diagonal with units on diagonal
                    D_inv_diag = []
                    for i in range(snf_basis.r):
                        d_ii = snf_basis._D[i, i]
                        D_inv_diag.append(domain.exquo(domain.one, d_ii))
                    D_inv = Matrix.diag(*D_inv_diag)
                    inv = snf_basis._V * D_inv * snf_basis._U
            elif combined_basis.cols < combined_basis.rows:
                # Tall matrix (more rows than columns): has right inverse but no left inverse
                # This shouldn't happen in a well-formed module basis computation
                # Use left pseudoinverse: (A^T A)^{-1} A^T (only works over fields)
                if _is_field(domain):
                    ATA = combined_basis.T * combined_basis
                    ATA_inv = ATA.inv()
                    inv = ATA_inv * combined_basis.T
                else:
                    raise ValueError(
                        f"Cannot invert tall matrix over PID: "
                        f"{combined_basis.rows}x{combined_basis.cols}. "
                        f"Module basis should have full rank equal to ambient dimension."
                    )
            # No else clause - cols > rows already handled above

        return tuple(pivots_per_block + [inv])

    def __repr__(self):
        return f"SNFMatrix {self.tolist()}"


class AlgebraicModule:
    """
    Finitely generated module over a domain.

    Represents M = span(span_set) / relations(zero_set).

    For ZZ (PIDs): May have torsion elements
    For QQ (fields): Always free (no torsion)
    """

    verify = True

    def __init__(self, domain: Domain, dim: int, span_set: List[Matrix], zero_set: List[Matrix]):
        """
        Initialize module.

        Args:
            domain: Ring to compute over (ZZ, QQ, etc.)
            dim: Dimension of ambient space
            span_set: List of column vectors generating the module
            zero_set: List of column vectors representing relations (elements that equal zero)
        """
        self.domain, self.dim = domain, dim
        self.span_set, self.zero_set = span_set, zero_set
        self._span_matrix = SNFMatrix.from_columns(domain, dim, span_set)
        self._zero_matrix = SNFMatrix.from_columns(domain, dim, zero_set)
        if AlgebraicModule.verify:
            self._verify_zero_in_span()

    def _verify_zero_in_span(self):
        """Ensure each zero_set vector lies in span(span_set) - required for well-defined quotient."""
        for k in self.zero_set:
            if not self._span_matrix.column_spans(k):
                raise AssertionError("zero_set element not in span(span_set)")

    def classify(self, v: Matrix) -> str:
        """Return in_zero_span, in_span_not_zero, or outside_span for vector v."""
        assert v.cols == 1 and v.rows == self.dim

        if len(self.zero_set) == 0:
            # With no explicit zero-generators we interpret K = {0}.
            in_zero_span = all(self.domain.convert(v[i, 0]) == self.domain.zero for i in range(self.dim))
        else:
            in_zero_span = self._zero_matrix.column_spans(v)

        if len(self.span_set) == 0:
            # Likewise, an empty span_set defines the zero submodule.
            in_span = all(self.domain.convert(v[i, 0]) == self.domain.zero for i in range(self.dim))
        else:
            in_span = self._span_matrix.column_spans(v)

        if in_zero_span:
            return 'in_zero_span'
        if in_span:
            return 'in_span_not_zero'
        return 'outside_span'


class Differential:
    """
    Linear map between modules, defined on generators plus zero relations.

    Represents a homomorphism d: M → M' where M and M' are AlgebraicModules.
    Can be constructed with partial information (non_zero_info) and completed.
    """

    verify = True

    def __init__(self, from_module: AlgebraicModule, to_module: AlgebraicModule, non_zero_info: Optional[Dict[Matrix, Matrix]] = None):
        """
        Initialize differential.

        Args:
            from_module: Source module
            to_module: Target module
            non_zero_info: Dict mapping generators to their images (optional)
        """
        self.from_module = from_module
        self.to_module = to_module
        self.matrix = None
        self.non_zero_info = non_zero_info or {}

        if Differential.verify:
            assert self.from_module.domain == self.to_module.domain
        self.domain = self.from_module.domain

        # Filter out entries whose source is in K; optionally verify target in K'.
        filtered: dict = {}
        for k, v in self.non_zero_info.items():
            if Differential.verify:
                assert v.rows == self.to_module.dim and v.cols == 1
                assert self.to_module._span_matrix.column_spans(v)
            if self.from_module._zero_matrix.column_spans(k):
                if Differential.verify:
                    assert self.to_module._zero_matrix.column_spans(v)
                continue
            filtered[k] = v
        self.non_zero_info = filtered

    def _solve_coeffs(self, basis: SNFMatrix, target: Matrix) -> Matrix:
        """Find one coefficient vector c with basis * c = target."""
        if basis.cols == 0:
            return Matrix.zeros(0, 1)
        if _is_field(self.domain):
            sol, _ = basis.gauss_jordan_solve(target)
            return sol
        # PID path via SNF
        U, D, V = basis._U, basis._D, basis._V
        r = basis.r
        bprime = U * target
        y = Matrix.zeros(basis.cols, 1)
        for i in range(r):
            di = D[i, i]
            y[i, 0] = self.domain.exquo(bprime[i, 0], di)
        x = V * y
        return x

    def _basis_and_redundancy(self, vectors: list[Matrix]):
        """Return independent subset, their indices, and representations for redundant ones."""
        basis_vecs: list[Matrix] = []
        basis_indices: list[int] = []
        reprs: dict[int, Matrix] = {}
        basis_mat = SNFMatrix.from_columns(self.domain, self.from_module.dim, [])

        for idx, v in enumerate(vectors):
            if basis_mat.cols > 0 and basis_mat.column_spans(v):
                coeffs = self._solve_coeffs(basis_mat, v)
                reprs[idx] = coeffs
                continue
            basis_vecs.append(v)
            basis_indices.append(idx)
            basis_mat = SNFMatrix.from_columns(self.domain, self.from_module.dim, basis_vecs)
        return basis_vecs, basis_indices, reprs

    def compute_matrix(self) -> bool:
        """
        Build the differential matrix if generators suffice.

        Returns False if non_zero_info together with K does not span S.
        On success, sets self.matrix (m' x m) and prunes redundant generators.
        """
        zero_set = list(self.from_module.zero_set)
        key_list = list(self.non_zero_info.keys())
        combined = zero_set + key_list

        # Do the keys ∪ K span S?
        generator_mat = SNFMatrix.from_columns(self.domain, self.from_module.dim, combined)
        for v in self.from_module.span_set:
            if not generator_mat.column_spans(v):
                return False

        # Extract a basis and record representations of redundant generators.
        basis_vecs, basis_indices, reprs = self._basis_and_redundancy(combined)
        basis_mat = SNFMatrix.from_columns(self.domain, self.from_module.dim, basis_vecs)

        # Build images for basis vectors.
        images_basis: list[Matrix] = []
        for idx in basis_indices:
            if idx < len(zero_set):
                images_basis.append(Matrix.zeros(self.to_module.dim, 1))
            else:
                key = key_list[idx - len(zero_set)]
                images_basis.append(self.non_zero_info[key])

        # Verify redundant images if requested.
        if Differential.verify:
            for idx, coeffs in reprs.items():
                img = Matrix.zeros(self.to_module.dim, 1)
                for j, c in enumerate(coeffs):
                    if c != 0:
                        img += images_basis[j] * c
                if idx < len(zero_set):
                    expected = Matrix.zeros(self.to_module.dim, 1)
                else:
                    key = key_list[idx - len(zero_set)]
                    expected = self.non_zero_info[key]
                assert img == expected

        # Compute matrix columns so that M * basis_vec = images_basis.
        m_rows = self.to_module.dim
        m_cols = self.from_module.dim
        result = Matrix.zeros(m_rows, m_cols)
        for j in range(m_cols):
            e = Matrix.zeros(m_cols, 1)
            e[j, 0] = self.domain.one if hasattr(self.domain, "one") else 1
            if basis_mat.column_spans(e):
                coeffs = self._solve_coeffs(basis_mat, e)
                col_img = Matrix.zeros(m_rows, 1)
                for k, c in enumerate(coeffs):
                    if c != 0:
                        col_img += images_basis[k] * c
                result[:, j] = col_img
            else:
                # Outside S; map to zero
                result[:, j] = Matrix.zeros(m_rows, 1)

        # Prune generators to the basis
        self.from_module.zero_set = [zero_set[i] for i in basis_indices if i < len(zero_set)]
        new_info = {}
        for idx in basis_indices:
            if idx >= len(zero_set):
                key = key_list[idx - len(zero_set)]
                new_info[key] = self.non_zero_info[key]
        self.non_zero_info = new_info

        self.matrix = result
        return True


class SpectralPageModule(AlgebraicModule):
    """Specialized module for spectral sequence pages E_r = ker(d_r) / im(d_{r-1}).

    Explicitly represents the quotient structure of a spectral sequence page:
    - Full module: Kernel of outgoing differential d_r
    - Relations: Image of incoming differential d_{r-1}

    This makes the quotient structure explicit for spectral sequence computations.

    Attributes:
        incoming_image: Basis of im(d_{r-1}) (becomes zero_genset)
        outgoing_kernel: Basis of ker(d_r) (defines ambient module)
    """

    def __init__(
        self,
        domain: Domain,
        dim: int,
        incoming_image: Optional[List[Matrix]] = None,
        outgoing_kernel: Optional[List[Matrix]] = None,
        page: Optional["Page"] = None,
        bidegree: Optional["Bidegree"] = None,
    ):
        """Initialize spectral page module E_r = ker(d_r) / im(d_{r-1}).

        Args:
            domain: ZZ or QQ
            dim: Dimension of ambient free module
            incoming_image: List of vectors forming basis of im(d_{r-1})
            outgoing_kernel: List of vectors forming basis of ker(d_r)
            page: Optional Page object containing this module (for context)
            bidegree: Optional Bidegree of this module

        Mathematical interpretation:
        - incoming_image are the relations (elements that die)
        - outgoing_kernel defines the full module before quotienting
        - The quotient ker/im gives the next page E_{r+1}
        """
        incoming_image = incoming_image or []
        outgoing_kernel = outgoing_kernel or []

        # Initialize base module with incoming image as relations
        super().__init__(domain, dim, span_set=outgoing_kernel, zero_set=incoming_image)

        # Store kernel explicitly for spectral sequence context
        self.outgoing_kernel = outgoing_kernel
        self._sp_basis_cache: Optional[Matrix] = None
        self._ker_basis_cache: Optional[Matrix] = None
        self.page = page
        self.bidegree = bidegree

    @property
    def sp_basis(self) -> Matrix:
        """Basis of non-trivial elements (surviving in quotient).

        Returns: Matrix whose columns are basis vectors of ker(d_r)/im(d_{r-1}).

        This represents the "true" elements of E_r that aren't killed by relations.
        """
        if self._sp_basis_cache is not None:
            return self._sp_basis_cache

        # Compute columnspace of spanning set to get basis
        if not self.span_set:
            self._sp_basis_cache = Matrix.zeros(self.dim, 0)
            return self._sp_basis_cache

        span_matrix = SNFMatrix.from_columns(self.domain, self.dim, self.span_set)
        self._sp_basis_cache = span_matrix.get_columnspace()
        return self._sp_basis_cache

    @property
    def ker_basis(self) -> Matrix:
        """Basis of killed elements (in im(d_{r-1})).

        Returns: Matrix whose columns are basis vectors of im(d_{r-1}).

        These elements are trivial in E_r (they die from previous differential).
        """
        if self._ker_basis_cache is not None:
            return self._ker_basis_cache

        if not self.zero_set:
            self._ker_basis_cache = Matrix.zeros(self.dim, 0)
            return self._ker_basis_cache

        # Compute columnspace of zero relations
        zero_matrix = SNFMatrix.from_columns(self.domain, self.dim, self.zero_set)
        self._ker_basis_cache = zero_matrix.get_columnspace()
        return self._ker_basis_cache

    @property
    def basis(self) -> Matrix:
        """Combined basis of kernel and spanning parts.

        Returns: Matrix whose columns are ker_basis followed by sp_basis.

        This is used by differential computation to express results in the full module basis.
        """
        return self.ker_basis.row_join(self.sp_basis)

    @property
    def abs_dim(self) -> int:
        """Absolute dimension (total dimension of ambient space).

        Returns: Dimension of the full ambient free module before quotienting.
        """
        return self.dim

    def classify_element(self, v: Matrix) -> Tuple[Matrix, Matrix]:
        """Classify element v in the spectral page module.

        Wrapper around inherited classify() for clarity in spectral sequence context.

        Args:
            v: Vector to classify (column vector of dimension self.dim)

        Returns: (span_coords, zero_coords) where:
        - span_coords: Coordinates of v in non-trivial basis (sp_basis)
        - zero_coords: Coordinates of v in killed basis (ker_basis)

        Interpretation:
        - If zero_coords != 0: v is trivial in E_r (killed by relations)
        - If span_coords != 0 and zero_coords = 0: v is non-trivial in E_r
        """
        return self.classify(v)




class Module:
    """Bigraded module for spectral sequence pages.

    Represents a finitely-generated module M = span(sp_basis) / zero_set where:
    - sp_basis: Basis for non-trivial elements (spanning basis)
    - zero_set (ker_basis): Elements that are identified with zero (kernel relations)

    This is the concrete module class used by Page objects for spectral sequence computation.
    Uses SNF machinery for exact arithmetic over both fields and PIDs.
    """

    def __init__(self, page: Page, bidegree: Bidegree, basis: Matrix, ker_basis: Matrix):
        """Initialize module with spanning and kernel bases.

        Args:
            page: Page object containing this module
            bidegree: Bidegree of this module in the spectral sequence
            basis: Matrix whose columns are spanning basis vectors (absolute coordinate space)
            ker_basis: Matrix whose columns are kernel basis vectors (relations)
        """
        self.page = page
        self.bidegree = bidegree

        # Use SNF-based multi_reduction to separate ker_basis and sp_basis
        ker_basis_idx, sp_basis_idx, self.basis_inv = SNFMatrix.multi_reduction_pid(
            page.ss.domain, ker_basis, basis
        )

        self.sp_basis: Matrix = basis[:, sp_basis_idx] if sp_basis_idx else Matrix.zeros(basis.rows, 0)
        self.ker_basis: Matrix = ker_basis[:, ker_basis_idx] if ker_basis_idx else Matrix.zeros(ker_basis.rows, 0)

    @property
    def basis(self):
        """Combined basis: ker_basis columns followed by sp_basis columns."""
        return self.ker_basis.row_join(self.sp_basis)

    @property
    def abs_dim(self):
        """Absolute dimension of ambient coordinate space."""
        return self.page.ss.get_abs_dimension(self.bidegree)

    @property
    def dim(self):
        """Dimension of the quotient module (number of non-trivial generators)."""
        return self.sp_basis.cols

    def __contains__(self, e: HomoElem):
        """Check if element e is in this module."""
        if e.page != self.page:
            return False
        if e.isZero():
            return True
        return e.bidegree == self.bidegree

    def classify(self, vec) -> int:
        """Classify a coordinate vector in this module.

        Args:
            vec: Column vector in absolute coordinate space

        Returns:
            0 if vec is in kernel (relations)
            1 if vec is in module but not kernel (non-trivial)
            2 if vec is outside module (error condition)
        """
        if self.abs_dim == 0:
            if vec.is_zero_matrix:
                return 0
            else:
                return 2

        assert vec.cols == 1 and vec.rows == self.abs_dim

        # Transform to basis coordinates
        indicator = self.basis_inv * vec

        # Find the last non-zero entry
        i = len(indicator) - 1
        domain_zero = self.page.ss.domain.zero if hasattr(self.page.ss.domain, 'zero') else self.page.ss.domain(0)
        while i >= 0 and self.page.ss.domain.convert(indicator[i, 0]) == domain_zero:
            i -= 1

        # Classify based on position in basis
        ker_size = self.ker_basis.cols
        sp_size = self.sp_basis.cols

        if i >= ker_size + sp_size:
            return 2  # Not spanned by basis
        elif i >= ker_size:
            return 1  # In sp_basis but not ker_basis
        else:
            return 0  # In ker_basis


class Differential:
    """Differential map between modules in a spectral sequence.

    Represents a map d: M → M' where M and M' are bigraded modules.
    Stores partial knowledge (non_zero_info) about where generators map,
    then derives the matrix representation using exact linear algebra.

    Uses SNF machinery to handle both field and PID computations correctly.
    """

    def __init__(self, page: Page, io_pairs: dict, d_bigrade: Bidegree):
        """Initialize differential with HomoElem-based knowledge.

        Args:
            page: Page object containing this differential
            io_pairs: Dict mapping element expressions to their images
            d_bigrade: Bidegree shift of the differential

        The io_pairs dict maps generator expressions to their images.
        This is stored as HomoElem pairs internally.
        """
        from src.element import HomoElem

        self.page = page
        self.d_bidegree = d_bigrade

        # Store HomoElem knowledge for compatibility with element.py
        self.knowledge: dict = {}
        for key, value in io_pairs.items():
            self.knowledge[HomoElem(page, key)] = HomoElem(page, value)

        # Cache for computed matrices
        self.calculated_matrices = SortedDict()

    def get_matrix(self, bidegree: Bidegree, debug=False):
        """Compute differential matrix at given bidegree using SNF machinery.

        This implements the full Leibniz rule-based differential computation,
        using SNF for exact arithmetic over both fields and PIDs.

        Args:
            bidegree: Source bidegree for the differential
            debug: If True, print debug information

        Returns:
            Matrix representation of differential at bidegree
        """
        if bidegree in self.calculated_matrices:
            return self.calculated_matrices[bidegree]

        from src.element import HomoElem, Bidegree as BidegreeType

        target_bidegree = bidegree + self.d_bidegree
        if self.page.ss.get_abs_dimension(target_bidegree) == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bidegree)])

        target_module = self.page[target_bidegree]
        target_dim = target_module.dim
        if target_dim == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bidegree)])

        # Get all known elements with correct bidegree
        pre_basis = Matrix()
        keys_to_combine = []
        bideg_to_combine = []
        temp_knowledge = {}

        for e in self.knowledge.keys():
            if e.bidegree == bidegree:
                pre_basis = pre_basis.row_join(e.coordinate)
                temp_knowledge[e] = self.knowledge[e]

                # Verify target bidegree is correct
                assert self.knowledge[e].bidegree is None or \
                       self.knowledge[e].bidegree == e.bidegree + self.d_bidegree

            elif not e.isZero():
                keys_to_combine.append(e)
                bideg_to_combine.append(e.bidegree)

        # Apply Leibniz rule to generate new differentials from known ones
        if len(bideg_to_combine) > 0 and bidegree != BidegreeType([0, 0]):
            combos = convex_integral_combinations(Matrix.hstack(*bideg_to_combine), bidegree)
            for combo in combos:
                if debug:
                    print(keys_to_combine, bideg_to_combine, combo, "target", bidegree)
                element = HomoElem(self.page, expr="1")
                d_value = HomoElem(self.page, expr="0")
                for i, exponent in enumerate(combo):
                    cur_key = keys_to_combine[i]
                    element *= cur_key ** exponent

                    # Leibniz rule: d(xy) = x*d(y) + y*d(x)
                    # Here applied as: d(x^n) = n*x^(n-1)*d(x)
                    temp = HomoElem(self.page, exponent) * cur_key ** (exponent - 1) * self.knowledge[cur_key]
                    for j, _exponent in enumerate(combo):
                        if j != i:
                            temp *= keys_to_combine[j] ** _exponent
                    d_value += temp
                if element.isZero():
                    continue

                if debug:
                    print(f"Generated Differential: d({element}) = {d_value}")

                temp_knowledge[element] = d_value
                pre_basis = pre_basis.row_join(element.coordinate)

                assert d_value.bidegree is None or \
                       d_value.bidegree == element.bidegree + self.d_bidegree, \
                       str(d_value.bidegree) + str(element.bidegree) + str(self.d_bidegree)

        # Expand pre_basis to a full basis of the source module using SNF machinery
        module = self.page.get_module(bidegree)
        target_module = self.page[target_bidegree]

        # Use SNF-aware multi_reduction for exact computation over domain
        # We need to span the quotient module, which has basis = ker_basis | sp_basis
        pivots_pre, pivots_module, inv = SNFMatrix.multi_reduction_pid(
            self.page.ss.domain, pre_basis, module.basis
        )

        res = Matrix([[]] * target_dim)

        # Process columns from pre_basis (known differentials)
        for i in pivots_pre:
            elem = HomoElem(self.page, abs_bideg=bidegree, abs_coordinate=pre_basis.col(i))
            target = temp_knowledge[elem]
            if target.isZero():
                a = Matrix([0] * target_dim)
                res = res.row_join(a)
            else:
                res = res.row_join(target.coordinate)

        # Process columns from module.basis (need to compute or ask for differential)
        # These are elements that are part of the module basis but not yet known
        for i in pivots_module:
            # Get the actual column from module.basis
            basis_col = module.basis.col(i)

            # Check if this column is in ker_basis or sp_basis
            if i < module.ker_basis.cols:
                # This is a kernel element (from im(d_{n-1}) at previous page)
                # For the differential, we must check if it's zero or ask user
                unknown_elem = HomoElem(self.page, abs_bideg=bidegree, abs_coordinate=basis_col)
                if debug:
                    print(f"Page {self.page.page_num}:\n"
                          f"\tUnknown differential for kernel element at {bidegree.tolist()}")

                if unknown_elem.isZero():
                    a = Matrix([0] * target_dim)
                    res = res.row_join(a)
                else:
                    expr = input(f"Please input d_{self.page.page_num}({unknown_elem})")
                    target = HomoElem(self.page, expr=expr)
                    if target.isZero():
                        a = Matrix([0] * target_dim)
                        res = res.row_join(a)
                    else:
                        res = res.row_join(target.coordinate)
            else:
                # This is from sp_basis (the non-trivial spanning part)
                unknown_idx = i - module.ker_basis.cols
                unknown_elem = HomoElem(self.page, abs_bideg=bidegree, abs_coordinate=module.sp_basis.col(unknown_idx))

                if debug:
                    print(f"Page {self.page.page_num}:\n"
                          f"\tOne unknown differential at abs_coordinate {module.sp_basis.col(unknown_idx).tolist()} "
                          f"encountered in the spanning basis of module {bidegree.tolist()}. The absolute basis here is "
                          f"{self.page.ss.get_abs_basis(bidegree)}")

                # If element is zero, don't ask for its differential
                if unknown_elem.isZero():
                    a = Matrix([0] * target_dim)
                    res = res.row_join(a)
                else:
                    expr = input(f"Please input d_{self.page.page_num}({unknown_elem})")
                    target = HomoElem(self.page, expr=expr)
                    if target.isZero():
                        a = Matrix([0] * target_dim)
                        res = res.row_join(a)
                    else:
                        res = res.row_join(target.coordinate)

        # Transform the differential matrix to coordinates in the target module basis
        # res has columns of d(basis_elem_i) in absolute target coordinates
        # We need to express these in target_module's basis coordinates
        output = target_module.basis_inv * res
        self.calculated_matrices[bidegree] = output
        return output

    def __call__(self, e):
        """Apply differential to element (stub)."""
        pass
