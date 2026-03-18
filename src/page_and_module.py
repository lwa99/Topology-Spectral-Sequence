from __future__ import annotations
from typing import TYPE_CHECKING

from src.snf import *
from src.matrices import *
from src.differential import Differential
from src.element import Bidegree, HomoElem, HomoCollection
from collections.abc import Iterable

if TYPE_CHECKING:
    from src.spectral_sequence import SpectralSequence

_verify = True


class Module:
    def __init__(self, page: Page, bidegree, span_set: Iterable[DMatrix], relation_set: Iterable[DMatrix]):
        self.page = page
        self.bideg = bidegree
        self.domain = self.page.domain
        self.dim = None  # TODO: compute dim from bideg
        self.span = HomoCollection(page=page, bideg=bidegree, coords=span_set)
        self.S = self.span.to_SNF_matrix()
        self.relation = HomoCollection(page=page, bideg=bidegree, coords=relation_set)
        self.R = self.relation.to_SNF_matrix()
        print(f"module initialization: bidegree:{bidegree}, span_set: {span_set}, S: {self.S}, R: {self.R}")

        if _verify:
            if self.S is None:
                for r in self.relation.coords:
                    # If the span set is empty, only literal zero vectors are valid relations.
                    assert all(e == self.domain.zero for e in r.to_list_flat())
            elif self.S.shape[0] == 0:
                for r in self.relation.coords:
                    # In a zero ambient module, only zero-dimensional relation vectors are valid.
                    assert r.shape[0] == 0
            else:
                for r in self.relation.coords:
                    assert self.S.spans(r)

    def get_structural_information(self):
        """
        Return a list of generators and their corresponding torsion information.
        Free part is followed by torsion part and marked by torsion 0.
        """
        if self.S is None:
            return None

        if self.R is None:
            return self.S.columns(), [self.domain.zero] * self.S.shape[1]

        P, Q, D = SNF.align(self.R, self.S)
        gens = (self.S * Q).columns()
        diag = D.diagonal()
        torsion = [diag[i] if i < len(diag) else self.domain.zero for i in range(len(gens))]
        return gens, torsion

    def classify(self, v: DMatrix):
        """Classify a coordinate vector in this module.

        Args:
            v: Column vector in absolute coordinate space

        Returns:
            0 if vec is in zero space (relations)
            1 if vec is in module but not zero space (non-trivial)
            2 if vec is outside module (error condition)
        """
        for i in v.to_list():
            if i != self.domain.zero:
                break
        else:
            return 0
        if self.S is None or not self.S.spans(v):
            return 2
        if self.R is None or not self.R.spans(v):
            return 1
        return 0

    def get_diff_span(self, I: HomoCollection = None, d_I: HomoCollection = None):
        """
        Return d(S) for this module's span generators.

        The computation now lives in Differential.complete_info_set()/get_diff_span(),
        so I and d_I are ignored and kept only for backward compatibility.
        """
        return self.page.d.get_diff_span(self.bideg)

    def get_diff_ker(self):
        """
        Compute ker(d) in absolute coordinates at this bidegree, then include source relations.

        If d(S) is represented by columns and target relations are R', we solve
            d(S) * u - R' * v = 0
        by taking the kernel of [d(S) | -R'] and projecting to the u-part.
        """
        # If the source span is empty, the kernel is exactly the relation part.
        if self.span.is_empty:
            return self.relation

        target_bideg = self.bideg + self.page.d.d_bidegree
        target_module = self.page[target_bideg]
        dS = self.page.d.get_diff_span(self.bideg)
        dS_M = dS.to_matrix()

        # Number of source span generators (columns of S).
        n = self.S.shape[1]

        if dS_M is None:
            # This can only occur for the trivial source span; guarded above, but keep safe.
            ker_coeff = DMatrix.from_list([[] for _ in range(n)], self.domain)
        else:
            if target_module.R is None:
                block = dS_M
            else:
                block = DMatrix.static_hstack(dS_M, -target_module.R)

            ker_block = SNF.kernel_of(block)
            ker_coeff = ker_block.extract(list(range(n)), list(range(ker_block.shape[1])))

        ker_from_span = self.S * ker_coeff
        ker_collection = HomoCollection.from_matrix(self.page, self.bideg, ker_from_span)
        return ker_collection.join(self.relation)


class Page:
    """
    Spectral sequence page.
    SpectralSequence.add_page()` constructs this object.
    """

    ss: "SpectralSequence"
    page_num: int

    def __init__(self, ss: "SpectralSequence", page_num: int, io_pairs: dict, d_bigrade: "Bidegree"):
        self.ss: SpectralSequence = ss
        self.domain = self.ss.domain
        self.page_num = page_num
        self.modules: dict = {}
        self.d = Differential(self, io_pairs, Bidegree(d_bigrade))

    @staticmethod
    def _normalize_bidegree(bidegree) -> Bidegree:
        if isinstance(bidegree, IMatrix):
            if bidegree.shape == (2, 1):
                return bidegree
            if bidegree.shape == (1, 2):
                return Bidegree([bidegree[0, 0], bidegree[0, 1]])
            raise ValueError(f"Bidegree matrix must have shape (2, 1) or (1, 2), got {bidegree.shape}.")

        if isinstance(bidegree, tuple) and len(bidegree) == 2:
            return Bidegree(bidegree)
        if isinstance(bidegree, list) and len(bidegree) == 2:
            return Bidegree(bidegree)
        raise TypeError(f"Invalid bidegree type: {type(bidegree)}. Expected 2-vector or 2-entry tuple/list.")

    def __getitem__(self, bidegree: Bidegree) -> Module:
        bidegree = self._normalize_bidegree(bidegree)
        if bidegree in self.modules:
            return self.modules[bidegree]
        else:
            output = self.generate_module(bidegree)
            self.modules[bidegree] = output
            return output

    def generate_module(self, bidegree) -> Module:
        if self.page_num == 1:
            d = self.ss.get_abs_dimension(bidegree)
            identity = DMatrix.eye((d, d), self.domain).columns()
            relations = self.ss.get_ker_basis(bidegree)
            return Module(self, bidegree, identity, relations)

        prev_page = self.ss.pages[self.page_num - 1]
        assert prev_page is not None

        # Ambient module for E_{r+1} is ker(d_r) inside E_r at this bidegree.
        prev_module_at_bideg = prev_page[bidegree]
        outgoing_kernel = prev_module_at_bideg.get_diff_ker()

        # Relations for E_{r+1}: image of incoming d_r from bidegree - d_r.
        incoming_source_bideg = bidegree - prev_page.d.d_bidegree
        if self.ss.get_abs_dimension(incoming_source_bideg) == 0:
            incoming_image = HomoCollection(page=prev_page, bideg=bidegree, coords=[])
        else:
            incoming_image = prev_page.d.get_diff_span(incoming_source_bideg)
            incoming_target_bideg = incoming_image.bideg
            assert incoming_target_bideg == bidegree, (
                f"Incoming differential target mismatch while building page {self.page_num} "
                f"at bidegree {bidegree}: expected target {bidegree}, got {incoming_target_bideg} "
                f"from page {prev_page.page_num}."
            )

        # Keep previous-page relations as zero in absolute coordinates of the new page.
        relations = incoming_image.join(prev_module_at_bideg.relation)
        return Module(self, bidegree, outgoing_kernel.coords, relations.coords)

    def divide(self, x: HomoElem, y: HomoElem):
        """
        Find q such that xq = y in the target module.

        Args:
            x: divider element
            y: dividend element

        Returns:
            (q_or_none, K_or_none) where:
            - q_or_none is one particular solution (None if unsolvable)
            - K_or_none generates the quotient ambiguity span(K_raw)/span(R_q),
              where K_raw is the full solution-difference kernel in absolute
              coordinates and R_q are relations in the quotient module.
        """
        q_bideg = y.bidegree - x.bidegree
        M_q, M_y = self[q_bideg], self[y.bidegree]
        xS_q = x * M_q.span
        xS_q_with_rel = xS_q.join(M_y.relation)

        if xS_q_with_rel.is_empty:
            if y.isZero():
                q_abs_dim = self.ss.get_abs_dimension(q_bideg)
                q_zero_coord = DMatrix.zeros((q_abs_dim, 1), self.domain)
                q = HomoElem(self, abs_bideg=q_bideg, abs_coordinate=q_zero_coord)
                K = DMatrix.zeros((q_abs_dim, 0), self.domain)
                return q, K
            else:
                return None, None

        solve_res = SNF.solve(y.coordinate, xS_q_with_rel.to_matrix())
        if solve_res is None:
            return None, None
        combined_coord, ker = solve_res

        source_col_num = len(xS_q)
        if source_col_num == 0:
            q_abs_dim = self.ss.get_abs_dimension(q_bideg)
            q_zero_coord = DMatrix.zeros((q_abs_dim, 1), self.domain)
            q = HomoElem(self, abs_bideg=q_bideg, abs_coordinate=q_zero_coord)
        else:
            coord_in_S = combined_coord.to_list()[:source_col_num]
            abs_coord = M_q.S * DMatrix.from_list(coord_in_S, self.domain)
            q = HomoElem(self, abs_bideg=q_bideg, abs_coordinate=abs_coord)

        if source_col_num == 0:
            q_abs_dim = self.ss.get_abs_dimension(q_bideg)
            K_raw = DMatrix.zeros((q_abs_dim, 0), self.domain)
        else:
            ker_coeff = ker.extract(list(range(source_col_num)), list(range(ker.shape[1])))
            K_raw = M_q.S * ker_coeff
        K_mod_rel = self._kernel_mod_relations(K_raw, M_q)
        return q, K_mod_rel

    def _kernel_mod_relations(self, K_raw: DMatrix, M_q: Module) -> DMatrix:
        """
        Reduce raw kernel ambiguity by quotient-module relations.

        Input:
            K_raw: generators of solution differences in absolute coordinates.
            M_q:   module where the quotient lives.

        Return:
            Generators representing span(K_raw)/span(R_q), where R_q = M_q.R.
        """
        if K_raw.shape[1] == 0:
            return K_raw
        if M_q.R is None or M_q.R.shape[1] == 0:
            return K_raw

        # Expected in well-defined module actions: quotient relations are
        # among solution-difference directions.
        solved = SNF.solve(M_q.R, K_raw)  # M_q.R = K_raw * X
        assert solved is not None, (
            "Expected quotient-module relations to lie in the divide-kernel ambiguity, "
            f"but failed for bidegree {M_q.bideg} on page {self.page_num}."
        )
        X = solved[0]
        _, Q, D = SNF.align(M_q.R, K_raw, _X=X)
        K_aligned = K_raw * Q

        diag_len = min(D.shape)
        keep_indices = []
        for i in range(K_aligned.shape[1]):
            if i < diag_len and self.domain.is_unit(D[i, i].element):
                # Unit diagonal means this direction is already in R_q.
                continue
            keep_indices.append(i)

        if len(keep_indices) == 0:
            return DMatrix.zeros((K_raw.shape[0], 0), self.domain)
        return K_aligned.extract_columns(keep_indices)

    def collection_divide_by(self, X: HomoCollection, l: list[HomoElem]):
        """Divide the i-th element in X by the i-th element in l and form a new HomoCollection."""
        elems = [self.divide(l[i], x)[0] for (i, x) in enumerate(X.elems)]
        return HomoCollection(page=self, bideg=X.bideg, elems=elems)
