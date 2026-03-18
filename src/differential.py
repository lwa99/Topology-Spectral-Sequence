from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING
from sympy.polys.polyerrors import ExactQuotientFailed

from src.snf import *
from src.matrices import *
from src.element import HomoElem, HomoCollection, Bidegree

if TYPE_CHECKING:
    from src.page_and_module import Page


@dataclass(frozen=True)
class DiffDivisionEdge:
    """Edge metadata for divisibility among known differential sources."""
    parent: HomoElem
    child: HomoElem
    quotient: HomoElem
    is_unique: bool


class DiffInfo:
    """
    Differential knowledge base.

    Stores known source->target differential values and maintains a divisibility
    graph over known sources. If a divides b, then a is recorded as a child of b.

    Notes:
    - This class is add-only.
    - Divisibility is currently checked by exact polynomial quotient.
    """

    def __init__(self, differential: Differential):
        self.differential = differential
        self.page = differential.page
        self.domain = differential.domain
        self._info: dict[HomoElem, HomoElem] = {}
        self._info_by_bideg: dict[Bidegree, dict[HomoElem, HomoElem]] = {}
        self._children: dict[HomoElem, dict[HomoElem, DiffDivisionEdge]] = {}
        self._parents: dict[HomoElem, dict[HomoElem, DiffDivisionEdge]] = {}

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

    def children_of(self, parent: HomoElem):
        return list(self._children.get(parent, {}).values())

    def parents_of(self, child: HomoElem):
        return list(self._parents.get(child, {}).values())

    def edges(self):
        out = []
        for child_map in self._children.values():
            out.extend(child_map.values())
        return out

    def add(self, src: HomoElem, tgt: HomoElem) -> bool:
        """
        Add one known differential value and update the divisibility graph.

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
        self._update_division_graph_for(src)
        return True

    def _record_edge(self, parent: HomoElem, child: HomoElem, quotient: HomoElem, is_unique: bool):
        edge = DiffDivisionEdge(parent=parent, child=child, quotient=quotient, is_unique=is_unique)
        if parent not in self._children:
            self._children[parent] = {}
        self._children[parent][child] = edge
        if child not in self._parents:
            self._parents[child] = {}
        self._parents[child][parent] = edge

    def _exact_division(self, divisor: HomoElem, dividend: HomoElem):
        """Return (quotient, unique) if divisor divides dividend exactly, else None."""
        if divisor.isZero():
            return None
        if dividend.isZero():
            return HomoElem(self.page, 0), True
        try:
            q_poly = dividend.poly.exquo(divisor.poly)
        except ExactQuotientFailed:
            return None
        except Exception:
            return None
        return HomoElem(self.page, q_poly.as_expr()), True

    def _update_division_graph_for(self, new_src: HomoElem):
        known = list(self._info.keys())
        for other in known:
            if other is new_src:
                continue

            # new_src divides other => new_src is a child of other
            res = self._exact_division(new_src, other)
            if res is not None:
                q, unique = res
                self._record_edge(parent=other, child=new_src, quotient=q, is_unique=unique)

            # other divides new_src => other is a child of new_src
            res = self._exact_division(other, new_src)
            if res is not None:
                q, unique = res
                self._record_edge(parent=new_src, child=other, quotient=q, is_unique=unique)


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
        # Leibniz propagation can create cross-bidegree dependencies; clear all cached slices.
        self.info_collections.clear()
        self.diff_span_cache.clear()
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
        if bidegree in self.diff_span_cache:
            return True
        self.extend_by_forward_leibniz(bidegree)
        dS, missing = self._infer_diff_span_from_current_info(bidegree)
        if dS is None:
            return False
        self.diff_span_cache[bidegree] = dS
        return len(missing) == 0

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
            added_pairs.extend(self.extend_by_forward_leibniz(bidegree))
            dS, missing_sources = self._infer_diff_span_from_current_info(bidegree)
            if dS is not None:
                self.diff_span_cache[bidegree] = dS
                return added_pairs

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
