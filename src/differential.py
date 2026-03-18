from __future__ import annotations
from typing import TYPE_CHECKING

from src.snf import *
from src.matrices import *
from src.element import HomoElem, HomoCollection, Bidegree

if TYPE_CHECKING:
    from src.page_and_module import Page


class DiffInfo:
    """
    Differential knowledge base with a distinguished lowest-level source set.

    This class is add-only. It stores known source->target differential values,
    and keeps `lowest`: nontrivial known sources that are not currently recovered
    as products of known sources via divide/reverse-Leibniz inference.
    """

    def __init__(self, differential: Differential):
        self.differential = differential
        self.page = differential.page
        self.domain = differential.domain
        self._info: dict[HomoElem, HomoElem] = {}
        self._info_by_bideg: dict[Bidegree, dict[HomoElem, HomoElem]] = {}
        self._lowest: set[HomoElem] = set()
        self._reason: dict[HomoElem, str] = {}

    def __len__(self):
        return len(self._info)

    def __contains__(self, src: HomoElem):
        return src in self._info

    def get(self, src: HomoElem, default=None):
        return self._info.get(src, default)

    def items(self):
        return self._info.items()

    def items_at(self, bidegree: Bidegree):
        return self._info_by_bideg.get(bidegree, {}).items()

    def sources(self):
        return list(self._info.keys())

    def sources_at(self, bidegree: Bidegree):
        return list(self._info_by_bideg.get(bidegree, {}).keys())

    def targets_at(self, bidegree: Bidegree):
        return list(self._info_by_bideg.get(bidegree, {}).values())

    def lowest_sources(self):
        return list(self._lowest)

    def lowest_sources_at(self, bidegree: Bidegree):
        return [s for s in self._lowest if s.bidegree == bidegree]

    def reason_of(self, src: HomoElem):
        return self._reason.get(src)

    @staticmethod
    def _is_unique_kernel(K: DMatrix | None, domain) -> bool:
        if K is None:
            return False
        if K.shape[1] == 0:
            return True
        return all(x == domain.zero for x in K.to_list_flat())

    def _is_nontrivial(self, e: HomoElem) -> bool:
        if e.isZero():
            return False
        if e.poly.total_degree() != 0:
            return True
        terms = e.poly.terms()
        if len(terms) != 1:
            return False
        coef = terms[0][1]
        return not self.domain.is_unit(coef)

    def _kernel_gens(self, bidegree: Bidegree, K: DMatrix | None):
        if K is None or K.shape[1] == 0:
            return []
        return HomoCollection.from_matrix(self.page, bidegree, K).elems

    def _constant_content_nonunit(self, e: HomoElem):
        """
        Return a non-unit common coefficient of all terms of e, if any.

        For example over ZZ:
          2*a -> 2
          6*a + 4*b -> 2
          a + 2*b -> None
        """
        if e.isZero():
            return None
        terms = e.poly.terms()
        if len(terms) == 0:
            return None
        c = terms[0][1]
        for _, coef in terms[1:]:
            c = self.domain.gcd(c, coef)
        if self.domain.is_ZZ and c < 0:
            c = -c
        if c == self.domain.zero or self.domain.is_unit(c):
            return None
        return c

    def _infer_from_constant_coeff(self, src: HomoElem, tgt: HomoElem, _seen: set[HomoElem]):
        """
        If src = c * q with non-unit scalar c and both divides are unique,
        infer d(q) = tgt / c and cache it recursively.
        """
        c = self._constant_content_nonunit(src)
        if c is None:
            return

        c_elem = HomoElem(self.page, self.domain.to_sympy(c))
        q_src, K_src = self.page.divide(c_elem, src)
        if q_src is None or not self._is_unique_kernel(K_src, self.domain):
            return
        if tgt.isZero():
            q_tgt = HomoElem(self.page, 0)
        else:
            q_tgt, K_tgt = self.page.divide(c_elem, tgt)
            if q_tgt is None or not self._is_unique_kernel(K_tgt, self.domain):
                return

        self.add_known(
            q_src,
            q_tgt,
            reason=f"const_divide:{c}",
            propagate=True,
            _seen=_seen,
        )

    def _solve_unique(self, left: HomoElem, rhs: HomoElem, *, rhs_bideg: Bidegree | None = None):
        if rhs.isZero() and rhs_bideg is not None:
            abs_dim = self.page.ss.get_abs_dimension(rhs_bideg)
            zero_coord = DMatrix.zeros((abs_dim, 1), self.domain)
            rhs_use = HomoElem(self.page, abs_bideg=rhs_bideg, abs_coordinate=zero_coord)
        else:
            rhs_use = rhs

        q, K = self.page.divide(left, rhs_use)
        if q is None:
            return None
        if not self._is_unique_kernel(K, self.domain):
            return None
        return q

    def _equal_up_to_relations(self, lhs: HomoElem, rhs: HomoElem, bidegree: Bidegree):
        diff = lhs - rhs
        if diff.isZero():
            return True
        return self.page[bidegree].classify(diff.coordinate) == 0

    def _recover_factors_via_left(self, product: HomoElem, left: HomoElem):
        """
        Try to recover d(q0) and d(k_j) from a known pair product = left * q.

        Return:
            (q0, dq0, kernel_pairs) on success, else None
            where kernel_pairs = [(k_j, d(k_j)), ...]
        """
        d_product = self._info.get(product)
        d_left = self._info.get(left)
        if d_product is None or d_left is None:
            return None

        q0, Kq = self.page.divide(left, product)
        if q0 is None:
            return None
        if not self._is_nontrivial(left) or not self._is_nontrivial(q0):
            return None

        target_bideg = product.bidegree + self.differential.d_bidegree
        rhs0 = d_product - (d_left * q0)
        dq0 = self._solve_unique(left, rhs0, rhs_bideg=target_bideg)
        if dq0 is None:
            return None

        q_bideg = product.bidegree - left.bidegree
        kernel_pairs: list[tuple[HomoElem, HomoElem]] = []
        for k in self._kernel_gens(q_bideg, Kq):
            if k.isZero():
                continue
            rhs_k = HomoElem(self.page, 0) - (d_left * k)
            dk = self._solve_unique(left, rhs_k, rhs_bideg=target_bideg)
            if dk is None:
                return None
            kernel_pairs.append((k, dk))

        return q0, dq0, kernel_pairs

    def _try_factor_and_cache(self, product: HomoElem, left: HomoElem, _seen: set[HomoElem]):
        recovered = self._recover_factors_via_left(product, left)
        if recovered is None:
            return False
        q0, dq0, kernel_pairs = recovered

        self.add_known(
            q0,
            dq0,
            reason=f"reverse_leibniz:{left}|{product}",
            propagate=True,
            _seen=_seen,
        )
        for k, dk in kernel_pairs:
            self.add_known(
                k,
                dk,
                reason=f"reverse_leibniz_kernel:{left}|{product}",
                propagate=True,
                _seen=_seen,
            )
        return True

    def _update_from_new_known(self, src: HomoElem, _seen: set[HomoElem]):
        if not self._is_nontrivial(src):
            self._lowest.discard(src)
            return

        represented_by_known = False
        current_lowest = list(self._lowest)
        for left in current_lowest:
            if left is src:
                continue
            if left not in self._info:
                continue
            if self._try_factor_and_cache(src, left, _seen):
                represented_by_known = True

        if represented_by_known:
            self._lowest.discard(src)
        else:
            self._lowest.add(src)

        # Keep lowest as an anti-chain: if src can represent an old lowest
        # element through known factors, that old element is no longer lowest.
        for old in list(self._lowest):
            if old is src:
                continue
            if old not in self._info:
                self._lowest.discard(old)
                continue
            if self._try_factor_and_cache(old, src, _seen):
                self._lowest.discard(old)

    def add(self, src: HomoElem, tgt: HomoElem) -> bool:
        return self.add_known(src, tgt, reason="manual", propagate=True)

    def add_known(
        self,
        src: HomoElem,
        tgt: HomoElem,
        *,
        reason: str = "manual",
        propagate: bool = True,
        _seen: set[HomoElem] | None = None,
    ) -> bool:
        """
        Add one known differential value and optionally propagate inference.

        Returns:
            True iff this source is newly added.
        """
        self.differential._validate_io_pair(src, tgt)
        if src in self._info:
            if self._info[src] != tgt:
                raise ValueError(
                    f"Conflicting differential data on page {self.page.page_num} at source bidegree {src.bidegree} "
                    f"for source {src}: {self._info[src]} vs {tgt}."
                )
            return False

        self._info[src] = tgt
        if src.bidegree not in self._info_by_bideg:
            self._info_by_bideg[src.bidegree] = {}
        self._info_by_bideg[src.bidegree][src] = tgt
        self._reason[src] = reason

        if propagate:
            if _seen is None:
                _seen = set()
            if src in _seen:
                return True
            _seen.add(src)
            try:
                self._infer_from_constant_coeff(src, tgt, _seen)
                self._update_from_new_known(src, _seen)
            finally:
                _seen.remove(src)
        return True

    def query_d(self, src: HomoElem, _seen: set[HomoElem] | None = None):
        """
        Try to compute d(src) from known values by product decomposition
        against lowest-level elements.

        On success, caches the result through add_known(..., reason="query").
        """
        known = self._info.get(src)
        if known is not None:
            return known
        if src.isZero():
            return HomoElem(self.page, 0)

        if _seen is None:
            _seen = set()
        if src in _seen:
            return None

        _seen.add(src)
        try:
            target_bideg = src.bidegree + self.differential.d_bidegree
            for left in list(self._lowest):
                d_left = self._info.get(left)
                if d_left is None:
                    continue
                if not self._is_nontrivial(left):
                    continue

                q0, Kq = self.page.divide(left, src)
                if q0 is None:
                    continue
                if not self._is_nontrivial(q0):
                    continue

                dq0 = self.query_d(q0, _seen)
                if dq0 is None:
                    continue

                q_bideg = src.bidegree - left.bidegree
                dk_map: dict[HomoElem, HomoElem] = {}
                ok = True
                for k in self._kernel_gens(q_bideg, Kq):
                    dk = self.query_d(k, _seen)
                    if dk is None:
                        ok = False
                        break
                    dk_map[k] = dk
                if not ok:
                    continue

                candidate = (left * dq0) + (d_left * q0)
                for k, dk in dk_map.items():
                    alt = (left * (dq0 + dk)) + (d_left * (q0 + k))
                    if not self._equal_up_to_relations(candidate, alt, target_bideg):
                        ok = False
                        break
                if not ok:
                    continue

                self.add_known(src, candidate, reason="query", propagate=True, _seen=_seen)
                return self._info[src]
            return None
        finally:
            _seen.remove(src)


class Differential:
    def __init__(self, page: Page, io_pairs: dict, d_bigrade: Bidegree):
        self.page = page
        self.domain = self.page.domain
        self.d_bidegree = d_bigrade
        self.diff_info = DiffInfo(self)
        # Total knowledge base on this page:
        # src -> d(src), plus bidegree-indexed views for fast lookups.
        self.info: dict[HomoElem, HomoElem] = {}
        self.info_by_bideg: dict[Bidegree, dict[HomoElem, HomoElem]] = {}
        self.info_collections: dict[Bidegree, HomoCollection] = {}
        # Cache: source bidegree -> d(source span generators) as a HomoCollection.
        self.diff_span_cache: dict[Bidegree, HomoCollection] = {}

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
        # Keep mirrors synchronized in case DiffInfo inferred values recursively.
        self._merge_from_diff_info()
        if src in self.info:
            if self.info[src] != tgt:
                raise ValueError(
                    f"Conflicting differential data on page {self.page.page_num} at source bidegree {src.bidegree} "
                    f"for source {src}: {self.info[src]} vs {tgt}."
                )
            return False

        added = self.diff_info.add_known(src, tgt, reason="manual", propagate=True)
        self._merge_from_diff_info()
        return added

    def _merge_from_diff_info(self) -> bool:
        """
        Merge knowledge inferred/stored in DiffInfo into Differential mirrors.

        Returns:
            True iff any new source was merged.
        """
        changed = False
        for src, tgt in self.diff_info.items():
            if src in self.info:
                if self.info[src] != tgt:
                    raise ValueError(
                        f"Conflicting differential data on page {self.page.page_num} at source bidegree {src.bidegree} "
                        f"for source {src}: {self.info[src]} vs {tgt}."
                    )
                continue
            self.info[src] = tgt
            if src.bidegree not in self.info_by_bideg:
                self.info_by_bideg[src.bidegree] = {}
            self.info_by_bideg[src.bidegree][src] = tgt
            changed = True

        if changed:
            # New knowledge can affect all cached slices and inferred spans.
            self.info_collections.clear()
            self.diff_span_cache.clear()
        return changed

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

    @staticmethod
    def _dedupe_elems(elems: list[HomoElem]) -> list[HomoElem]:
        res = []
        seen = set()
        for e in elems:
            if e in seen:
                continue
            seen.add(e)
            res.append(e)
        return res

    @staticmethod
    def _kernel_unique_up_to_relations(K: DMatrix | None, module) -> bool:
        if K is None:
            return False
        for j in range(K.shape[1]):
            col = K.extract(list(range(K.shape[0])), [j])
            if module.classify(col) != 0:
                return False
        return True

    def _known_target_or_zero(self, src: HomoElem, zero_target: HomoElem, target_bideg: Bidegree):
        if src.isZero():
            return zero_target
        known = self.info.get(src)
        if known is None:
            return None
        if known.isZero():
            return zero_target
        assert known.bidegree == target_bideg, (
            f"Known d-image has unexpected bidegree on page {self.page.page_num} "
            f"for source {src}: got {known.bidegree}, expected {target_bideg}."
        )
        return known

    def _infer_diff_span_from_current_info(self, bidegree):
        """
        Try to infer d(S) at one source bidegree from current known info.

        Returns:
            (dS, missing_sources)
            - dS is HomoCollection if fully inferred, else None
            - missing_sources are transformed source generators requiring user input
        """
        module = self.page[bidegree]
        target_bideg = bidegree + self.d_bidegree
        target_module = self.page[target_bideg]

        if module.span.is_empty:
            return HomoCollection(page=self.page, bideg=target_bideg, coords=[]), []

        target_abs_dim = self.page.ss.get_abs_dimension(target_bideg)
        zero_coord = DMatrix.zeros((target_abs_dim, 1), self.domain)
        zero_target = HomoElem(self.page, abs_bideg=target_bideg, abs_coordinate=zero_coord)

        I, d_I, _ = self.info_with_images_at(bidegree)

        # No equations yet: ask/derive on the original source generators directly.
        if I.is_empty:
            missing = []
            y_elems: list[HomoElem | None] = [None] * module.S.shape[1]
            for i, src in enumerate(module.span.elems):
                known = self._known_target_or_zero(src, zero_target, target_bideg)
                if known is None:
                    if module.classify(src.coordinate) == 0:
                        y_elems[i] = zero_target
                    else:
                        missing.append(src)
                else:
                    y_elems[i] = known
            if missing:
                return None, self._dedupe_elems(missing)
            y = HomoCollection(page=self.page, bideg=target_bideg, elems=[e for e in y_elems if e is not None])
            return y, []

        I_M = I.to_matrix()
        P, Q, D = SNF.align(I_M, module.S)
        rhs = d_I * P
        transformed_sources = HomoCollection.from_matrix(self.page, bidegree, module.S * Q)

        n, m = D.shape
        diag_len = min(n, m)
        rhs_cols = rhs.elems

        y_elems: list[HomoElem | None] = [None] * n
        missing_sources: list[HomoElem] = []

        for i in range(diag_len):
            src_i = transformed_sources.elems[i]
            known = self._known_target_or_zero(src_i, zero_target, target_bideg)
            if known is not None:
                y_elems[i] = known
                continue

            d_i = D[i, i].element
            rhs_i = rhs_cols[i]
            if d_i == self.domain.zero:
                if not rhs_i.isZero():
                    raise ValueError(
                        f"Inconsistent differential alignment on page {self.page.page_num} "
                        f"at source bidegree {bidegree}: expected zero rhs at diagonal index {i}."
                    )
                if module.classify(src_i.coordinate) == 0:
                    y_elems[i] = zero_target
                else:
                    missing_sources.append(src_i)
                continue

            divider = HomoElem(self.page, self.domain.to_sympy(d_i))
            q, K = self.page.divide(divider, rhs_i)
            if q is None:
                raise ValueError(
                    f"Cannot divide inferred differential column by diagonal entry {d_i} "
                    f"on page {self.page.page_num} at source bidegree {bidegree}."
                )

            if self._kernel_unique_up_to_relations(K, target_module):
                y_elems[i] = q
            else:
                missing_sources.append(src_i)

        for j in range(diag_len, len(rhs_cols)):
            if not rhs_cols[j].isZero():
                raise ValueError(
                    f"Inconsistent differential alignment on page {self.page.page_num} "
                    f"at source bidegree {bidegree}: nonzero rhs column {j} "
                    f"outside diagonal part of D with shape {D.shape}."
                )

        # Rows beyond diagonal part are unconstrained by YD = rhs; they need direct values.
        for i in range(diag_len, n):
            src_i = transformed_sources.elems[i]
            known = self._known_target_or_zero(src_i, zero_target, target_bideg)
            if known is None:
                if module.classify(src_i.coordinate) == 0:
                    y_elems[i] = zero_target
                else:
                    missing_sources.append(src_i)
            else:
                y_elems[i] = known

        missing_sources = self._dedupe_elems(missing_sources)
        if missing_sources:
            return None, missing_sources

        for i in range(n):
            if y_elems[i] is None:
                y_elems[i] = zero_target

        y = HomoCollection(page=self.page, bideg=target_bideg, elems=[e for e in y_elems if e is not None])
        Q_inv = SNF.invert_unimodular(Q)
        return y * Q_inv, []

    def complete_info_set(self, bidegree):
        """
        Complete differential data at one source bidegree and cache d(S).

        It first tries to infer d(S) via the aligned equation:
            d(I) * P = d(S) * Q * D
        using division by diagonal entries of D. User input is only requested for
        transformed generators that remain ambiguous/unconstrained.
        """
        module = self.page[bidegree]
        target_bideg = bidegree + self.d_bidegree
        if module.span.is_empty:
            self.diff_span_cache[bidegree] = HomoCollection(page=self.page, bideg=target_bideg, coords=[])
            return []

        added_pairs: list[tuple[HomoElem, HomoElem]] = []
        while True:
            # Legacy forward-Leibniz expansion intentionally disabled:
            # added_pairs.extend(self.extend_by_forward_leibniz(bidegree))
            for src in module.span.elems:
                if src in self.info:
                    continue
                inferred = self.diff_info.query_d(src)
                if inferred is not None:
                    added_pairs.append((src, inferred))
            self._merge_from_diff_info()

            dS, missing_sources = self._infer_diff_span_from_current_info(bidegree)
            if dS is not None:
                self.diff_span_cache[bidegree] = dS
                return added_pairs

            inferred_any = False
            for src in missing_sources:
                if src in self.info:
                    continue
                if self.diff_info.query_d(src) is not None:
                    inferred_any = True
            if inferred_any:
                self._merge_from_diff_info()
                continue

            gens = ", ".join(str(g) for g in self.page.ss.gen)
            print(
                f"Need {len(missing_sources)} additional differential value(s) on page {self.page.page_num} "
                f"at bidegree {bidegree}. "
                f"Please enter each image using generators ({gens})."
            )

            new_added = 0
            for src in missing_sources:
                expr = input(
                    f"[page {self.page.page_num}, bidegree {bidegree}] Enter d({src}) [default 0]: "
                ).strip()
                if expr == "":
                    expr = "0"
                tgt = HomoElem(self.page, expr)
                if self._add_info_pair(src, tgt):
                    added_pairs.append((src, tgt))
                    new_added += 1

            if new_added == 0:
                raise ValueError(
                    f"Unable to complete differential data on page {self.page.page_num} at bidegree {bidegree}: "
                    f"no new information was added but d(S) is still ambiguous."
                )

    def get_diff_span(self, bidegree):
        if bidegree not in self.diff_span_cache:
            self.complete_info_set(bidegree)
        return self.diff_span_cache[bidegree]
