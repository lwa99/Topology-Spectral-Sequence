from sortedcontainers import SortedDict
from scalar import Scalar
from abstractPage import AbstractPage
from numpy import ndarray, array, asarray
from copy import copy, deepcopy


class DOArray(ndarray):
    """
    One dimensional array equipped with dictionary order. Used to store degrees of generators.
    """
    def __new__(cls, input_array):
        return asarray(input_array).view(cls)

    def __lt__(self, other):
        for i, n in enumerate(self):
            if n < other[i]:
                return True
        return False

    def __eq__(self, other):
        # for i, n in enumerate(self):
        #     if n != other[i]:
        #         return False
        # return True
        return all(super().__eq__(other))

    def __hash__(self):
        return int(sum(self))


class Element:
    """
    An element in a page is a polynomial of the generators.
    Coefficients in every polynomial are non-zero (zero terms would be deleted in real time).
    """
    def __init__(self, page: AbstractPage, degrees: ndarray = None, coef: Scalar = None):
        self.page = page
        self.coefMap = SortedDict()
        self.bigrade = None
        if degrees is not None:
            assert len(degrees) == self.page.gen_num
            self.coefMap[DOArray(degrees)] = coef
            self.bigrade = degrees * self.page.generator_bigrades

    def __add__(self, other: 'Element'):
        if len(self.coefMap) == 0:
            return deepcopy(other)
        output = deepcopy(self)
        if len(other.coefMap) == 0:
            return output

        for key, scalar in other.coefMap.items():

            # If the key is absent, the following line will add the key associated with an empty scalar wrapper.
            # Then, the value is replaced by the actual value.
            # IF the key is present, the following line will do nothing and return the current value, which will be
            # modified accordingly.
            # If it becomes 0, the key-value pair will be deleted.
            cur_val = output.coefMap.setdefault(key, self.page.get_scalar(None))
            if cur_val == self.page.get_scalar(None):
                cur_val.update(scalar.val)
            else:
                cur_val.increase_by(scalar)

                # Now take care of 0 coefficient
                if cur_val.val == 0:
                    del output.coefMap[key]

        # TODO: recalculate bigrade
        return output

    def __mul__(self, other):
        pass

    def _bigrade_from_deg_config(self, deg_config: DOArray | tuple[int, ...]):
        output = 0
        for i, n in enumerate(deg_config):
            output += n * self.page.generator_bigrades[i]
        return output

    def __str__(self):
        output = ""
        if len(self.coefMap) == 0:
            return "zero polynomial (element)"
        for key, value in self.coefMap.items():
            output += str(value.val)
            for i, degree in enumerate(key):
                output += f"({self.page.generators[i]}^{degree})"
            output += " + "
        output = output[:-2] + f" (mod {self.page.c})"
        return output


if __name__ == "__main__":
    dc_1 = DOArray((1, 2, 3))
    dc_2 = DOArray((1, 3, 2))
    for _i in dc_1:
        print(_i)
    print(dc_1 < dc_2)
