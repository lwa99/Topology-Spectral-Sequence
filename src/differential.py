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
        self.info_pairs: list[tuple[HomoElem, HomoElem]] = []
        self.info_collections: dict[Bidegree, HomoCollection] = {}

        # Safeguard against contradictory input-output data.
        for key, value in io_pairs.items():
            src = HomoElem(page, key)
            tgt = HomoElem(page, value)
            self._add_info_pair(src, tgt)
    
    def _validate_io_pair(self, src: HomoElem, tgt: HomoElem):
        """Validate one user-provided differential pair for obvious contradictions."""
        if src.isZero():
            if not tgt.isZero():
                raise ValueError(f"Invalid differential data: d(0) must be 0, but got {tgt}.")
            return

        target_bideg = src.bidegree + self.d_bidegree
        if not tgt.isZero() and tgt.bidegree != target_bideg:
            raise ValueError(
                "Invalid differential data bidegree mismatch: "
                f"d({src}) should land in bidegree {target_bideg}, but got {tgt.bidegree}."
            )

    def _add_info_pair(self, src: HomoElem, tgt: HomoElem) -> bool:
        """
        Insert one known differential value.

        Return:
            True iff a new source was added.
        """
        self._validate_io_pair(src, tgt)
        if src in self.info:
            if self.info[src] != tgt:
                raise ValueError(f"Conflicting differential data on source {src}: {self.info[src]} vs {tgt}.")
            return False

        self.info[src] = tgt
        self.info_pairs.append((src, tgt))
        self.info_collections.pop(src.bidegree, None)
        return True

    def info_collection_at(self, bidegree):
        """Gather the known information correspond to the specified bidegree."""
        if bidegree in self.info_collections.keys():
            return self.info_collections[bidegree]
        info_collection = HomoCollection(page=self.page, bideg=bidegree,
                                         coords=[src.coordinate for (src, _) in self.info_pairs
                                                 if src.bidegree == bidegree])
        res = info_collection.join(self.page[bidegree].relation)
        self.info_collections[bidegree] = res
        return res

    def info_with_images_at(self, bidegree):
        """
        Return I and d(I) at a bidegree, where I includes known info plus source relations.
        """
        I = self.info_collection_at(bidegree)
        target_bideg = bidegree + self.d_bidegree

        target_abs_dim = self.page.ss.get_abs_dimension(target_bideg)
        zero_coord = DMatrix.from_list([[self.domain.zero] for _ in range(target_abs_dim)], self.domain)
        zero_target = HomoElem(self.page, abs_bideg=target_bideg, abs_coordinate=zero_coord)

        dI_elems: list[HomoElem] = []
        for e in I.elems:
            image = None
            for src, tgt in self.info_pairs:
                if src.bidegree == e.bidegree and src.coordinate == e.coordinate:
                    image = tgt
                    break

            if image is None or image.isZero():
                dI_elems.append(zero_target)
            else:
                assert image.bidegree == target_bideg, \
                    f"Known d-image has unexpected bidegree: got {image.bidegree}, expected {target_bideg}"
                dI_elems.append(image)

        d_I = HomoCollection(page=self.page, bideg=target_bideg, elems=dI_elems)
        return I, d_I, target_bideg

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

    def complete_info_set(self, bidegree):
        """
        Expand known differential inputs at one bidegree minimally until they span the source module.

        If I = span(known inputs + source relations) and S = source span matrix, write I = S * X.
        Let U * X * V = D be SNF. Then I * V = (S * U^{-1}) * D.
        A diagonal entry d_i is a unit iff the i-th transformed source generator is already controlled by I.
        We only ask for d() on transformed generators with non-unit d_i.
        """
        module = self.page[bidegree]
        if module.span.is_empty:
            return []

        n = module.S.shape[1]
        I = self.info_collection_at(bidegree)

        if I.is_empty:
            missing_indices = list(range(n))
            transformed_sources = module.span
        else:
            I_M = I.to_matrix()
            solved = SNF.solve(I_M, module.S)
            if solved is None:
                raise ValueError(
                    f"Known information at {bidegree} does not lie in source span; cannot complete info set."
                )
            X = solved[0]
            D, U, _ = SNF.decomp(X)
            U_inv = U.inv_den()[0]
            transformed_sources = HomoCollection.from_matrix(self.page, bidegree, module.S * U_inv)

            diag_len = min(D.shape[0], D.shape[1])
            missing_indices = []
            for i in range(n):
                if i < diag_len:
                    d_i = D[i, i].element
                else:
                    d_i = self.domain.zero
                if not self.domain.is_unit(d_i):
                    missing_indices.append(i)

        if not missing_indices:
            return []

        gens = ", ".join(str(g) for g in self.page.ss.gen)
        print(
            f"Need {len(missing_indices)} additional differential value(s) at bidegree {bidegree}. "
            f"Please enter each image using generators ({gens})."
        )

        added_pairs: list[tuple[HomoElem, HomoElem]] = []
        for idx in missing_indices:
            src = transformed_sources.elems[idx]
            expr = input(f"Enter d({src}) [default 0]: ").strip()
            if expr == "":
                expr = "0"
            tgt = HomoElem(self.page, expr)
            if self._add_info_pair(src, tgt):
                added_pairs.append((src, tgt))

        self.info_collections.pop(bidegree, None)
        return added_pairs

    def refine_info_set(self, bidegree):
        return self.complete_info_set(bidegree)

    def expand_info_set(self, bidegree):
        return self.complete_info_set(bidegree)
