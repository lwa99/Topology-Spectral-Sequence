from __future__ import annotations
from typing import TYPE_CHECKING

from src.snf import *
from src.matrices import *
from src.element import HomoElem, HomoCollection, Bidegree

if TYPE_CHECKING:
    from src.page_and_module import Page


class Differential:
    def __init__(self, page: Page, io_pairs: dict, d_bigrade: Bidegree):
        self.page = page
        self.domain = self.page.domain
        self.d_bidegree = d_bigrade
        # Total knowledge base on this page:
        # src -> d(src), plus bidegree-indexed views for fast lookups.
        self.info: dict[HomoElem, HomoElem] = {}
        self.info_by_bideg: dict[Bidegree, dict[HomoElem, HomoElem]] = {}
        self.info_pairs: list[tuple[HomoElem, HomoElem]] = []
        self.info_collections: dict[Bidegree, HomoCollection] = {}

        # Safeguard against contradictory input-output data.
        for key, value in io_pairs.items():
            src = HomoElem(page, key)
            tgt = HomoElem(page, value)
            self._add_info_pair(src, tgt)

        # In a unital setting, d(1) = 0 should always hold.
        unit_src = HomoElem(page, 1)
        unit_tgt = HomoElem(page, 0)
        self._add_info_pair(unit_src, unit_tgt)

    def _validate_io_pair(self, src: HomoElem, tgt: HomoElem):
        """Validate one user-provided differential pair for obvious contradictions."""
        if src.isZero():
            if not tgt.isZero():
                raise ValueError(
                    f"Invalid differential data on page {self.page.page_num} at bidegree {src.bidegree}: "
                    f"d(0) must be 0, but got {tgt}."
                )
            return

        target_bideg = src.bidegree + self.d_bidegree
        if not tgt.isZero() and tgt.bidegree != target_bideg:
            raise ValueError(
                f"Invalid differential data bidegree mismatch on page {self.page.page_num} "
                f"at source bidegree {src.bidegree}: "
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
                raise ValueError(
                    f"Conflicting differential data on page {self.page.page_num} at source bidegree {src.bidegree} "
                    f"for source {src}: {self.info[src]} vs {tgt}."
                )
            return False

        self.info[src] = tgt
        if src.bidegree not in self.info_by_bideg:
            self.info_by_bideg[src.bidegree] = {}
        self.info_by_bideg[src.bidegree][src] = tgt
        self.info_pairs.append((src, tgt))
        # Leibniz propagation can create cross-bidegree dependencies; clear all cached slices.
        self.info_collections.clear()
        return True

    def info_collection_at(self, bidegree):
        """Gather the known information correspond to the specified bidegree."""
        if bidegree in self.info_collections.keys():
            return self.info_collections[bidegree]
        known_at_bideg = self.info_by_bideg.get(bidegree, {})
        info_collection = HomoCollection(page=self.page, bideg=bidegree,
                                         coords=[src.coordinate for src in known_at_bideg.keys()])
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
        zero_coord = DMatrix.zeros((target_abs_dim, 1), self.domain)
        zero_target = HomoElem(self.page, abs_bideg=target_bideg, abs_coordinate=zero_coord)

        dI_elems: list[HomoElem] = []
        for e in I.elems:
            image = self.info.get(e)

            if image is None or image.isZero():
                dI_elems.append(zero_target)
            else:
                assert image.bidegree == target_bideg, \
                    (
                        f"Known d-image has unexpected bidegree on page {self.page.page_num} "
                        f"at source bidegree {bidegree}: got {image.bidegree}, expected {target_bideg}"
                    )
                dI_elems.append(image)

        d_I = HomoCollection(page=self.page, bideg=target_bideg, elems=dI_elems)
        return I, d_I, target_bideg

    def factor_from_known(self, src: HomoElem, *, allow_unit_factor: bool = False):
        """
        Try to write src as a product of two known-source elements.

        Return:
            (left, right) if found, else None.
        """
        if src.isZero() or src.bidegree is None:
            return None

        for left_bideg, left_dict in self.info_by_bideg.items():
            right_bideg = src.bidegree - left_bideg
            right_dict = self.info_by_bideg.get(right_bideg)
            if right_dict is None:
                continue

            for left in left_dict.keys():
                if left.isZero():
                    continue
                if not allow_unit_factor and left.poly.total_degree() == 0:
                    continue

                for right in right_dict.keys():
                    if right.isZero():
                        continue
                    if not allow_unit_factor and right.poly.total_degree() == 0:
                        continue
                    if left * right == src:
                        return left, right
        return None

    def derive_by_leibniz(self, src: HomoElem):
        """If src factors through known sources, derive d(src) via Leibniz."""
        decomposition = self.factor_from_known(src)
        if decomposition is None:
            return None
        left, right = decomposition
        d_left = self.info.get(left)
        d_right = self.info.get(right)
        if d_left is None or d_right is None:
            return None
        return (left * d_right) + (d_left * right)

    def extend_by_forward_leibniz(self, bidegree):
        """
        Iteratively add new known d-values at one bidegree using Leibniz from existing knowledge.
        """
        module = self.page[bidegree]
        if module.span.is_empty:
            return []

        added: list[tuple[HomoElem, HomoElem]] = []
        changed = True
        while changed:
            changed = False
            for src in module.span.elems:
                if src in self.info:
                    continue
                derived = self.derive_by_leibniz(src)
                if derived is None:
                    continue
                if self._add_info_pair(src, derived):
                    added.append((src, derived))
                    changed = True
        return added

    def info_complete(self, bidegree):
        self.extend_by_forward_leibniz(bidegree)
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

        self.extend_by_forward_leibniz(bidegree)
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
                    f"Known information on page {self.page.page_num} at bidegree {bidegree} "
                    f"does not lie in source span; cannot complete info set."
                )
            X = solved[0]
            D, U, _ = SNF.decomp(X)
            U_inv = SNF.invert_unimodular(U)
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
            f"Need {len(missing_indices)} additional differential value(s) on page {self.page.page_num} "
            f"at bidegree {bidegree}. "
            f"Please enter each image using generators ({gens})."
        )

        added_pairs: list[tuple[HomoElem, HomoElem]] = []
        for idx in missing_indices:
            src = transformed_sources.elems[idx]
            expr = input(
                f"[page {self.page.page_num}, bidegree {bidegree}] Enter d({src}) [default 0]: "
            ).strip()
            if expr == "":
                expr = "0"
            tgt = HomoElem(self.page, expr)
            if self._add_info_pair(src, tgt):
                added_pairs.append((src, tgt))

        # User input may unlock further sources via Leibniz at the same bidegree.
        added_pairs.extend(self.extend_by_forward_leibniz(bidegree))
        self.info_collections.clear()
        return added_pairs

    def refine_info_set(self, bidegree):
        return self.complete_info_set(bidegree)

    def expand_info_set(self, bidegree):
        return self.complete_info_set(bidegree)
