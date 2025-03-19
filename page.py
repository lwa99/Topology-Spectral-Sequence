from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from module import Module
    from spectral_sequence import SpectralSequence

from utilities import Matrix, Vector
from sortedcontainers import SortedDict

# class Page(AbstractPage):
#     def set_differential(self):
#         """
#         Set the information of the differential
#         """
#
#         raise NotImplementedError
#
#     def get_scalar(self, val) -> Scalar:
#         return Scalar(self.c, val)
#
#     def get_module(self, bigrade: Matrix):
#         for module in self.calculated_modules.keys():
#             if bigrade == module:
#                 return module
#             return Module(self, bigrade)


class Page:
    def __init__(self, ss, page_num):
        self.ss: SpectralSequence = ss
        self.page_num = page_num
        self.subspaces = SortedDict([])

    def get_subspace(self, bigrade):
        if bigrade in self.subspaces:
            return self.subspaces[bigrade]
        else:
            output = self.generate_subspace(bigrade)
            self.subspaces[bigrade] = output
            return output

    def generate_subspace(self, bigrade) -> Module:
        from module import Module

        # This is the function that handles generates subspaces of new pages.
        if self.page_num == 1:
            return Module(self, bigrade,
                          Matrix.eye(self.ss.get_abs_dimension(bigrade)),
                          self.ss.get_ker_basis(bigrade)
                          )

        # TODO: build subspace from previous page
        raise NotImplementedError

    # def test_zero_homo_poly(self, coef_map: SortedDict):
    #     assert len(Vector) == len(self.ss.generators)
    #     first_exponent = coef_map.keys().__iter__().__next__()
    #     bigrade = self.ss.get_bigrade(first_exponent)
    #     subspace = self.get_subspace(bigrade)
    #     return subspace.kernelContains(self.ss.get_std_coordinate(coef_map))


if __name__ == "__main__":
    pass
