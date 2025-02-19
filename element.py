from sortedcontainers import SortedDict
from scalar import Scalar
from abstractPage import AbstractPage


class _DegreeConfiguration:
    """
    A degree configuration is an ordered list of degree on each generator. It is equipped with the dictionary
    order so that we can use it as a key in red black tree.
    """

    def __init__(self, degree_config: tuple[int]):
        self.degree_config = degree_config

    def __lt__(self, other):
        for i, n in enumerate(self.degree_config):
            if n < other.degree_config[i]:
                return True
        return False

    def __gt__(self, other):
        for i, n in enumerate(self.degree_config):
            if n > other.degree_config[i]:
                return True
        return False

    def __eq__(self, other):
        for i, n in enumerate(self.degree_config):
            if n != other.degree_config[i]:
                return False
        return True

    def __hash__(self):
        return sum(self.degree_config)


class Element:
    """
    An element in a page is a polynomial of the generators.
    """
    def __init__(self, page: AbstractPage, degrees=None, coef: Scalar =None):
        self.page = page
        self.coefMap = SortedDict()
        if degrees is not None:
            assert len(degrees) == self.page.gen_num
            self.coefMap[_DegreeConfiguration(degrees)] = coef

    def __add__(self, other: 'Element'):
        output = Element(self.page)
        output.coefMap = self.coefMap
        for key, scalar in other.coefMap.items():

            # If the key is absent, the following line will add the key associated with an empty scalar wrapper.
            # Then, the value is replaced by the actual value.
            # IF the key is present, the following line will do nothing and return the current value, which will be
            # modified accordingly.
            # If it becomes 0, the key-value pair will be deleted.
            cur_val = output.coefMap.setdefault(key, self.page.get_scalar(None))
            print(cur_val)
            if cur_val == self.page.get_scalar(None):
                cur_val.update(scalar.val)
            else:
                print("Add:" + str(scalar))
                cur_val.increase_by(scalar)
                if cur_val.val == 0:
                    del output.coefMap[key]
        return output

    def __str__(self):
        output = ""
        if len(self.coefMap) == 0:
            return "zero polynomial (element)"
        for key, value in self.coefMap.items():
            output += str(value.val)
            for i, degree in enumerate(key.degree_config):
                output += f"({self.page.generators[i]}^{degree})"
            output += " + "
        output = output[:-2] + f" (mod {self.page.c})"
        return output

    @property
    def bigrade(self):
        return
