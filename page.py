from scalar import Scalar
from abstractPage import AbstractPage
from element import Element
from module import Module
from utilities import Matrix, pprint


class Page(AbstractPage):
    def set_differential(self):
        """
        Set the information of the differential
        """

        raise NotImplementedError

    def get_scalar(self, val) -> Scalar:
        return Scalar(self.c, val)

    def get_module(self, bigrade: Matrix):
        for module in self.calculated_modules.keys():
            if bigrade == module:
                return module
            return Module(self, bigrade)

    def _calculate_next_page_module(self):
        """
        Calculate the module of the next page at a give bigrade
        """
        # Step 1: Locate the "previous" module and the "next" module



        # Step 2: Calculate the matrix representations of differential


        # Step 3: Return kernel / image
        raise NotImplementedError

    def _compare_elements(self, e1, e2):
        raise NotImplementedError


if __name__ == "__main__":
    try:
        assert False
    except AssertionError:
        print("Debug Mode: ON")

    p = Page(["x", "y", "z"], Matrix([[1, 2, 3], [3, 4, 5]]), 7, 1)
    print("Generator Bigrades: (Displayed as Columns)")
    pprint(p.generator_bigrades)
    e_1 = Element(p, Matrix([1, 2, 0]), p.get_scalar(3))
    print("e_1: ", e_1)
    e_2 = Element(p, Matrix([3, 1, 0]), p.get_scalar(5))
    e_3 = e_1 * e_2
    print("e_3: ", e_3)
    print("e_3 bigrade:")
    pprint(e_3.bigrade)
