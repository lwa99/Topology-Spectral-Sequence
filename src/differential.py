from __future__ import annotations
from typing import TYPE_CHECKING

from snf import *
from matrices import *
from element import HomoElem, HomoCollection, Bidegree

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

        self.info_collections: dict[Bidegree, HomoCollection] = {}

    def info_collection_at(self, bidegree):
        """Gather the known information correspond to the specified bidegree."""
        if bidegree in self.info_collections.keys():
            return self.info_collections[bidegree]
        info_collection = HomoCollection(page=self.page, bideg=bidegree,
                                         coords=[k.coordinate for k in self.info.keys() if k.bidegree == bidegree])
        res = info_collection.join(self.page[bidegree].relation)
        self.info_collections[bidegree] = res
        return res

    def extend_by_forward_leibniz(self, bidegree):
        pass

    def info_complete(self, bidegree):
        module = self.page[bidegree]
        if self.info_collection_at(bidegree).is_empty:
            return module.span.is_empty
        for s in module.span.coords:
            if not self.info_collection_at(bidegree).to_SNF_matrix().spans(s):
                return False
        return True

    def expand_info_set(self):
        pass
