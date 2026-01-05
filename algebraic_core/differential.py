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

    def I_at(self, bidegree):
        info_set: list[DMatrix] = [k.coordinate for k in self.info.keys() if k.bidegree == bidegree]
        info_set.extend(self.page[bidegree].relation_set)
        return SNFMatrix.hstack(*info_set, domain=self.domain)

    def extend_by_forward_leibniz(self, bidegree):
        pass

    def info_complete(self, bidegree):
        module = self.page[bidegree]
        K = SNFMatrix.hstack(self.I, module.relation_set, domain=self.domain)
        for s in module.span_set:
            if not K.spans(s):
                return False
        return True

    def expand_info_set(self):
        pass