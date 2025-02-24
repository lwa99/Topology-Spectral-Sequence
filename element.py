from sortedcontainers import SortedDict
from scalar import Scalar
from abstractPage import AbstractPage
from tools import Matrix
from copy import deepcopy


class Element:
    """
    An element in a page is a polynomial of the generators.
    Coefficients in every polynomial are non-zero (zero terms would be deleted in real time).

    Important Note:
    For zero polynomial, we define its bigrade to be empty: DOArray([]).
    For inhomogeneous (and thus non-zero) polynomial, we define its bigrade to be None.
    """
    def __init__(self, page: AbstractPage, degrees: Matrix = None, coef: Scalar = None):
        self.page = page
        self.coefMap = SortedDict()
        self.bigrade = Matrix([])
        if degrees is not None:
            assert len(degrees) == self.page.gen_num
            self.coefMap[degrees] = coef
            self.bigrade = self.page.generator_bigrades * degrees

    def __add__(self, other: 'Element'):
        if len(self.coefMap) == 0:
            return deepcopy(other)

        # From now on we may assume that self is not the zero polynomial
        output = deepcopy(self)
        if len(other.coefMap) == 0:
            return output
        term_deleted = False

        for key, scalar in other.coefMap.items():
            # If the key is absent, the following line will add the key associated with an empty scalar wrapper.
            # Then, the value is replaced by the actual value.
            # IF the key is present, the following line will do nothing and return the current value, which will be
            # modified accordingly.
            # If it becomes 0, the key-value pair will be deleted.
            cur_coef = output.coefMap.setdefault(key, self.page.get_scalar(None))
            if cur_coef == self.page.get_scalar(None):
                cur_coef.update(scalar.val)
            else:
                cur_coef.increase_by(scalar)

                # Now take care of 0 coefficient
                if cur_coef.val == 0:
                    del output.coefMap[key]
                    term_deleted = True

        # Now we update bigrades
        if output.bigrade is not None:
            # If the original bigrade exists but is not equal to that of the other, the bigrade become undefined
            if output.bigrade != other.bigrade:
                output.bigrade = None
        elif term_deleted:
            # If the original bigrade did not exist and adding a new polynomial results in deletion of some terms,
            # the resulting polynomial may become homogeneous again
            output._update_bigrade()

        return output

    def __mul__(self, other):
        if len(self.coefMap) == 0:
            return deepcopy(self)
        if len(other.coefMap) == 0:
            return other

        # From now on we may assume that both operands are non-zero
        output = Element(self.page)
        if self.bigrade is None or other.bigrade is None:
            output.bigrade = None
        else:
            output.bigrade = self.bigrade + other.bigrade

        for deg_1, coef_1 in self.coefMap.items():
            for deg_2, coef_2 in other.coefMap.items():
                output.coefMap[deg_1 + deg_2] = coef_1 * coef_2

        return output

    def _update_bigrade(self):
        if len(self.coefMap) == 0:
            self.bigrade = Matrix([])
            return
        _iter = self.coefMap.keys().__iter__()  # get the key iterator
        res = self.page.generator_bigrades * _iter.__next__()  # get the first key and calculate bigrade
        try:
            while True:
                if self.page.generator_bigrades * _iter.__next__() != res:
                    self.bigrade = None
                    return
        except StopIteration:
            assert isinstance(res, Matrix)
            self.bigrade = res
            return

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
    dc_1 = Matrix((1, 2, 3))
    dc_2 = Matrix((1, 3, 2))
    for _i in dc_1:
        print(_i)
    print(dc_1 < dc_2)
