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
            if self.S is None or self.S.shape[0] == 0:
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
            return self.S.columns, [self.domain.zero] * self.S.shape[1]

        P, Q, D = SNF.align(self.R, self.S)
        return (self.S * Q).columns, D.diagonal() + [self.domain.zero] * (self.S.shape[1] - self.R.shape[1])

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

    def get_diff_span(self, I: HomoCollection, d_I: HomoCollection):
        """
        Given I and d(I), compute d(S)
        """
        if not self.page.d.info_complete(self.bideg):
            self.page.d.complete_info_set(self.bideg)
            I, d_I, _ = self.page.d.info_with_images_at(self.bideg)
        if I.is_empty:
            return HomoCollection(coords=[], bideg=d_I.bideg, page=d_I.page)
        I_M = I.to_matrix()
        P, Q, D = SNF.align(I_M, self.S)
        # P, Q, D = self.S.align(I)
        print(D.diagonal())
        rhs = d_I * P

        # If D has shape (n, m), we solve Y * D = rhs where Y has n columns.
        # Columns beyond min(n,m) in rhs must be zero (D has zero columns there).
        n, m = D.shape
        diag_len = min(n, m)
        rhs_cols = rhs.elems
        target_abs_dim = self.page.ss.get_abs_dimension(rhs.bideg)
        zero_coord = DMatrix.zeros((target_abs_dim, 1), self.domain)
        zero_elem = HomoElem(self.page, abs_bideg=rhs.bideg, abs_coordinate=zero_coord)

        y_known_elems: list[HomoElem] = []
        for i in range(diag_len):
            d_i = D[i, i].element
            rhs_i = rhs_cols[i]
            if d_i == self.domain.zero:
                if not rhs_i.isZero():
                    raise ValueError(
                        f"Inconsistent differential alignment at {self.bideg}: expected zero rhs at diagonal index {i}."
                    )
                y_known_elems.append(zero_elem)
            else:
                divider = HomoElem(self.page, self.domain.to_sympy(d_i))
                q = self.page.divide(divider, rhs_i)
                if q is None:
                    raise ValueError(
                        f"Cannot divide differential image by diagonal entry at index {i} for bidegree {self.bideg}."
                    )
                y_known_elems.append(q)

        # Consistency for rhs columns not controlled by diagonal entries.
        for j in range(diag_len, len(rhs_cols)):
            if not rhs_cols[j].isZero():
                raise ValueError(
                    f"Inconsistent differential alignment at {self.bideg}: nonzero rhs column {j} "
                    f"outside diagonal part of D with shape {D.shape}."
                )

        # Extend Y with zero columns when D has more rows than columns.
        if n > diag_len:
            y_elems = y_known_elems + [zero_elem] * (n - diag_len)
            y = HomoCollection(page=self.page, bideg=rhs.bideg, elems=y_elems)
        else:
            y = HomoCollection(page=self.page, bideg=rhs.bideg, elems=y_known_elems)

        return y * Q.inv_den()[0]

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

        if not self.page.d.info_complete(self.bideg):
            self.page.d.complete_info_set(self.bideg)

        I, d_I, target_bideg = self.page.d.info_with_images_at(self.bideg)
        target_module = self.page[target_bideg]

        dS = self.get_diff_span(I, d_I)
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

    def __getitem__(self, bidegree: Bidegree) -> Module:
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
            incoming_source_module = prev_page[incoming_source_bideg]
            if not prev_page.d.info_complete(incoming_source_bideg):
                prev_page.d.complete_info_set(incoming_source_bideg)
            I_in, dI_in, incoming_target_bideg = prev_page.d.info_with_images_at(incoming_source_bideg)
            assert incoming_target_bideg == bidegree
            incoming_image = incoming_source_module.get_diff_span(I_in, dI_in)

        # Keep previous-page relations as zero in absolute coordinates of the new page.
        relations = incoming_image.join(prev_module_at_bideg.relation)
        return Module(self, bidegree, outgoing_kernel.coords, relations.coords)

    def divide(self, x: HomoElem, y: HomoElem, *, return_uniqueness: bool = False):
        """
        Find q such that xq = y in the target module.

        Args:
            x: divider element
            y: dividend element
            return_uniqueness: if True, also return whether q is unique up to relations

        Returns:
            - if return_uniqueness is False: q or None
            - if return_uniqueness is True: (q_or_none, is_unique_up_to_relations)
        """
        q_bideg = y.bidegree - x.bidegree
        M_q, M_y = self[q_bideg], self[y.bidegree]
        print("in divide", y, M_q.dim)
        xS_q = x * M_q.span
        xS_q_with_rel = xS_q.join(M_y.relation)

        if xS_q_with_rel.is_empty:
            if y.isZero():
                q_abs_dim = self.ss.get_abs_dimension(q_bideg)
                q_zero_coord = DMatrix.zeros((q_abs_dim, 1), self.domain)
                q = HomoElem(self, abs_bideg=q_bideg, abs_coordinate=q_zero_coord)
                if return_uniqueness:
                    return q, True
                return q
            else:
                if return_uniqueness:
                    return None, False
                return None

        solve_res = SNF.solve(y.coordinate, xS_q_with_rel.to_matrix())
        if solve_res is None:
            if return_uniqueness:
                return None, False
            return None
        combined_coord, ker = solve_res
        coord_in_S = combined_coord.to_list()[:len(xS_q)]
        abs_coord = M_q.S * DMatrix.from_list(coord_in_S, self.domain)
        q = HomoElem(self, abs_bideg=q_bideg, abs_coordinate=abs_coord)

        # Uniqueness up to relations:
        # q is unique mod relations iff every projected kernel direction is relation-trivial in M_q.
        is_unique_up_to_rel = True
        source_col_num = len(xS_q)
        for j in range(ker.shape[1]):
            delta_coeff = ker.extract(list(range(source_col_num)), [j])
            if source_col_num == 0:
                continue
            delta_abs = M_q.S * delta_coeff
            if M_q.classify(delta_abs) != 0:
                is_unique_up_to_rel = False
                break

        if return_uniqueness:
            return q, is_unique_up_to_rel
        return q

    def collection_divide_by(self, X: HomoCollection, l: list[HomoElem]):
        """Divide the i-th element in X by the i-th element in l and form a new HomoCollection."""
        elems = [self.divide(l[i], x) for (i, x) in enumerate(X.elems)]
        return HomoCollection(page=self, bideg=X.bideg, elems=elems)

