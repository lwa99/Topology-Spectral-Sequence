from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page

from utilities import Matrix, convex_integral_combinations
from element import HomoElem, Bidegree
from sortedcontainers import SortedDict
from sympy import Poly


class Differential:
    def __init__(self, page: Page, io_pairs: dict[Poly, Poly], d_bidegree: Bidegree):
        self.page = page
        self.d_bidegree = d_bidegree
        # knowledge maps domain element  -> image element
        self.knowledge: dict[HomoElem, HomoElem] = {
            HomoElem(page, k): HomoElem(page, v) for k, v in io_pairs.items()
        }
        # cache for already‑computed matrices, ordered by bidegree
        self.calculated_matrices: SortedDict[Bidegree, Matrix] = SortedDict()

    # ---------------------------------------------------------------------
    # Matrix of the differential on a fixed bidegree component
    # ---------------------------------------------------------------------
    def get_matrix(self, bigrade: Bidegree):
        """Return the matrix of *d* on the module E_r^{bigrade}. The result is
        cached so successive calls are cheap."""

        if bigrade in self.calculated_matrices:
            return self.calculated_matrices[bigrade]

        target_bidegree = bigrade + self.d_bidegree
        # If the target module is zero, the whole matrix is zero.
        if self.page.ss.get_abs_dimension(target_bidegree) == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bigrade)])

        target_dim = self.page[target_bidegree].dim
        if target_dim == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bigrade)])

        # -----------------------------------------------------------------
        # 1. Collect all *known* elements of the correct bidegree that we can
        #    use as a starting basis.
        # -----------------------------------------------------------------
        pre_basis = Matrix()
        keys_to_combine: list[HomoElem] = []
        bideg_to_combine: list[Bidegree] = []
        temp_knowledge: dict[HomoElem, HomoElem] = {}

        for e in self.knowledge:
            if e.bidegree == bigrade:
                pre_basis = pre_basis.row_join(e.coordinate)
                temp_knowledge[e] = self.knowledge[e]

                # image should either be zero or land in the correct shifted bidegree
                assert self.knowledge[e].bidegree is None or \
                       self.knowledge[e].bidegree == e.bidegree + self.d_bidegree
            elif not e.isZero():
                # keep for multiplicative generation later
                keys_to_combine.append(e)
                bideg_to_combine.append(e.bidegree)

        # -----------------------------------------------------------------
        # 2. Try to generate **new** relations multiplicatively when possible
        # -----------------------------------------------------------------
        if bideg_to_combine and bigrade != Bidegree([0, 0]):
            combos = convex_integral_combinations(Matrix.hstack(*bideg_to_combine), bigrade)
            for combo in combos:
                element = HomoElem(self.page, expr="1")
                d_value = HomoElem(self.page, expr="0")
                for i, exponent in enumerate(combo):
                    cur_key = keys_to_combine[i]
                    element *= cur_key ** exponent

                    temp = HomoElem(self.page, exponent) * cur_key ** (exponent - 1) * self.knowledge[cur_key]
                    for j, _exp in enumerate(combo):
                        if j != i:
                            temp *= keys_to_combine[j] ** _exp
                    d_value += temp
                if element.isZero():
                    continue

                # Record the generated differential in a *temporary* knowledge base.
                temp_knowledge[element] = d_value
                pre_basis = pre_basis.row_join(element.coordinate)

                assert d_value.bidegree is None or \
                       d_value.bidegree == element.bidegree + self.d_bidegree

        # -----------------------------------------------------------------
        # 3. Extend pre_basis to a full basis of the domain module and build
        #    the matrix column‑by‑column.
        # -----------------------------------------------------------------
        module = self.page.get_module(bigrade)
        pivots_1, pivots_2, pivots_3, inv = Matrix.multi_reduction(
            pre_basis, module.ker_basis, module.sp_basis
        )

        res = Matrix([[]] * target_dim)

        # (i) columns coming from known elements
        for i in pivots_1:
            elem = HomoElem(self.page, abs_bideg=bigrade, abs_coordinate=pre_basis.col(i))
            target = temp_knowledge[elem]
            res = res.row_join(Matrix.zeros(target_dim, 1) if target.isZero() else target.coordinate)

        # (ii) columns that were already in the kernel (zero image)
        for _ in pivots_2:
            res = res.row_join(Matrix.zeros(target_dim, 1))

        # (iii) columns whose image the program cannot deduce automatically
        for i in pivots_3:
            print(
                f"Page {self.page.page_num}: Unknown differential for basis vector {module.sp_basis.col(i).tolist()} in E_r^{bigrade}."
            )
            unknown_elem = HomoElem(self.page, abs_bideg=bigrade, abs_coordinate=module.sp_basis.col(i))
            expr = input(f"Please input d_{self.page.page_num}({unknown_elem}): ")
            target = HomoElem(self.page, expr=expr)
            res = res.row_join(Matrix.zeros(target_dim, 1) if target.isZero() else target.coordinate)

        # Convert to the standard basis of the domain module
        output = res * inv * module.basis
        self.calculated_matrices[bigrade] = output
        return output

    # ---------------------------------------------------------------------
    # Utility: insert new differential information
    # ---------------------------------------------------------------------
    def update_knowledge(self, iopair: dict[HomoElem, HomoElem]):
        for key, value in iopair.items():
            if key in self.knowledge and self.knowledge[key] != value:
                raise ValueError(
                    f"Conflict for key {key}: existing value {self.knowledge[key]}, new value {value}"
                )
            self.knowledge[key] = value
            # the image itself is always killed
            self.knowledge[value] = HomoElem(value.page, expr=0)

    # allow Differential(e) syntax
    def __call__(self, e: HomoElem):
        return self.knowledge.get(e, HomoElem(self.page, expr="0"))
