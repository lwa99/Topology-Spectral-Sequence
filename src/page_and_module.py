from __future__ import annotations
from typing import TYPE_CHECKING

from snf import *
from matrices import *
from differential import Differential
from element import Bidegree, HomoElem, HomoCollection
from collections.abc import Iterable

if TYPE_CHECKING:
    from spectral_sequence import SpectralSequence

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
        I_M = I.to_matrix()
        P, Q, D = SNF.align(I_M, self.S)
        # P, Q, D = self.S.align(I)
        print(D.diagonal())
        diag = [HomoElem(self.page, self.domain(x)) for x in D.diagonal()]
        assert len(diag) == D.shape[0] == D.shape[1]
        return self.page.collection_divide_by(d_I * P, diag) * Q.inv_den()[0]

    def get_next_module(self):
        pass


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

    def divide(self, x: HomoElem, y: HomoElem):
        """find q such that xq = y"""
        q_bideg = y.bidegree - x.bidegree
        M_q, M_y = self[q_bideg], self[y.bidegree]
        print("in divide", y, M_q.dim)
        xS_q = x * M_q.span
        xS_q_with_rel = xS_q.extend(M_y.relation)

        if xS_q_with_rel.is_empty:
            if y.isZero():
                return HomoElem(self, self.domain.zero)
            else:
                return None

        combined_coord = SNF.solve(y.coordinate, xS_q_with_rel.to_matrix())[0]
        if combined_coord is None:
            return None
        coord_in_S = combined_coord.to_list()[:len(xS_q)]
        abs_coord = M_q.S * DMatrix.from_list(coord_in_S, self.domain)
        return HomoElem(self, abs_bideg=q_bideg, abs_coordinate=abs_coord)

    def collection_divide_by(self, X: HomoCollection, l: list[HomoElem]):
        elems = [self.divide(l[i], x) for (i, x) in enumerate(X.elems)]
        return HomoCollection(elems=elems)

