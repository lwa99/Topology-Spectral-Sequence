from __future__ import annotations
from typing import TYPE_CHECKING

from snf import *
from matrices import *
from element import Bidegree, HomoElem

if TYPE_CHECKING:
    from page_and_module import Page


class Differential:
    def __init__(self, page: Page, io_pairs: dict, d_bigrade: Bidegree):
        self.page = page
        self.domain = self.page.domain
        self.d_bidegree = d_bigrade
        self.info: dict[HomoElem, HomoElem] = {}
        for key, value in io_pairs.items():
            self.info[HomoElem(page, key)] = HomoElem(page, value)

        self.info_matrices: dict[Bidegree, DMatrix] = {}

    def I_at(self, bidegree):
        """Gather the known information correspond to the specified bidegree."""
        info_set: list[DMatrix] = [k.coordinate for k in self.info.keys() if k.bidegree == bidegree]
        info_set.extend(self.page[bidegree].relation.elems)
        res = DMatrix.static_hstack(*info_set)
        self.info_matrices[bidegree] = res
        return res

    def extend_by_forward_leibniz(self, bidegree):
        pass

    def info_complete(self, bidegree):
        module = self.page[bidegree]
        K = SNFMatrix.static_hstack(self.info_matrices[bidegree], module.relation.to_matrix())
        for s in module.span.elems:
            if not K.spans(s):
                return False
        return True

    def expand_info_set(self):
        pass