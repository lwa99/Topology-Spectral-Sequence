from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from page import Page

from utilities import Matrix, convex_integral_combinations
from element import HomoElem, Bidegree
from sortedcontainers import SortedDict
from sympy import ZZ


class Differential:
    def __init__(self, page: Page, io_pairs: dict, d_bigrade: Bidegree):
        self.page = page
        self.d_bidegree = d_bigrade
        self.knowledge: dict[HomoElem, HomoElem] = {}
        for key, value in io_pairs.items():
            self.knowledge[HomoElem(page, key)] = HomoElem(page, value)
        self.calculated_matrices = SortedDict()

    def get_matrix(self, bidegree: Bidegree):
        if bidegree in self.calculated_matrices:
            return self.calculated_matrices[bidegree]

        target_bidegree = bidegree + self.d_bidegree
        if self.page.ss.get_abs_dimension(target_bidegree) == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bidegree)])

        target_dim = self.page[target_bidegree].dim
        if target_dim == 0:
            return Matrix([[0] * self.page.ss.get_abs_dimension(bidegree)])

        # Get all known elements with correct bidegree
        pre_basis = Matrix()
        keys_to_combine = []
        bideg_to_combine = []
        temp_knowledge = {}
        for e in self.knowledge.keys():
            if e.bidegree == bidegree:
                pre_basis = pre_basis.row_join(e.coordinate)
                temp_knowledge[e] = self.knowledge[e]

                # Make sure that the target bidegree is correct. It is either None (if the target is zero) or
                # it should be equal to the starting bidegree added by the differential bidegree.
                assert self.knowledge[e].bidegree is None or \
                       self.knowledge[e].bidegree == e.bidegree + self.d_bidegree

            elif not e.isZero():
                keys_to_combine.append(e)
                bideg_to_combine.append(e.bidegree)

        if len(bideg_to_combine) > 0 and bidegree != Bidegree([0, 0]):
            combos = convex_integral_combinations(Matrix.hstack(*bideg_to_combine), bidegree)
            # TODO: add a bound on exponents to utilize the torsion of elements (we don't care the case where the
            #  computed element here is zero.)
            #  Also put a cap of exponent of 2 on t if t^3 is present
            #  Also, maybe we should move this step
            for combo in combos:
                print(keys_to_combine, bideg_to_combine, combo, "target", bidegree)
                element = HomoElem(self.page, expr="1")
                d_value = HomoElem(self.page, expr="0")
                for i, exponent in enumerate(combo):
                    cur_key = keys_to_combine[i]
                    element *= cur_key ** exponent

                    temp = HomoElem(self.page, exponent) * cur_key ** (exponent - 1) * self.knowledge[cur_key]
                    for j, _exponent in enumerate(combo):
                        if j != i:
                            temp *= keys_to_combine[j] ** _exponent
                    d_value += temp
                if element.isZero():
                    continue

                print(f"Generated Differential: d({element}) = {d_value}")
                # self.knowledge[element] = d_value
                # This was deleted because we don't want the knowledge base to have multiplicative relations.

                temp_knowledge[element] = d_value
                pre_basis = pre_basis.row_join(element.coordinate)

                assert d_value.bidegree is None or \
                       d_value.bidegree == element.bidegree + self.d_bidegree, \
                       str(d_value.bidegree) + str(element.bidegree) + str(self.d_bidegree)


        '''
        In either cases: we start from the known elements, add the zero-space elements and then ask for new elements
        when necessary.
        '''
        if self.page.ss.domain == ZZ:
            '''
            In the Z case, we take 
            
            '''
        else:
            # Expand pre_basis to a basis of the starting module
            module = self.page.get_module(bidegree)
            pivots_1, pivots_2, pivots_3, inv = Matrix.multi_reduction(pre_basis, module.ker_basis, module.sp_basis)

            res = Matrix([[]] * target_dim)

            for i in pivots_1:
                elem = HomoElem(self.page, abs_bideg=bidegree, abs_coordinate=pre_basis.col(i))
                target = temp_knowledge[elem]
                if target.isZero():
                    a = Matrix([0] * target_dim)
                    res = res.row_join(a)
                else:
                    res = res.row_join(target.coordinate)

            for _ in pivots_2:
                # For the zero basis, their differential is zero, since they are already zero.
                res = res.row_join(Matrix.zeros(res.rows, 1))

            for i in pivots_3:
                print(f"Page {self.page.page_num}:\n"
                      f"\tOne unknown differential at abs_coordinate {module.sp_basis.col(i).tolist()} "
                      f"encountered in the spanning basis of module {bidegree.tolist()}. The absolute basis here is "
                      f"{self.page.ss.get_abs_genset(bidegree)}")
                unknown_elem = HomoElem(self.page, abs_bideg=bidegree, abs_coordinate=module.sp_basis.col(i))
                expr = input(f"Please input d_{self.page.page_num}({unknown_elem})")

                target = HomoElem(self.page, expr=expr)
                if target.isZero():
                    a = Matrix([0] * target_dim)
                    res = res.row_join(a)
                else:
                    res = res.row_join(target.coordinate)

            output = res * inv * module.basis
            self.calculated_matrices[bidegree] = output
            return output
