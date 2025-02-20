"""
Deprecated!!!!
Math Structures for the Project
"""
class BiGradedAlgebra:
    class DegreeConfiguration:
        """
        A degree configuration is an ordered list of degree on each generator. It is equipped with the dictionary
        order to speed up searching.
        """
        def __init__(self, degree_config: tuple[int]):
            self.degree_config = degree_config

        def __le__(self, other):
            for i, n in enumerate(self.degree_config):
                if n < other.degree_config[i]:
                    return True
            return False

        def __ge__(self, other):
            for i, n in enumerate(self.degree_config):
                if n > other.degree_config[i]:
                    return True
            return False

        def __eq__(self, other):
            for i, n in enumerate(self.degree_config):
                if n != other.degree_config[i]:
                    return False
            return True

    class Polynomial:
        """
        Homogeneous Term

        Every polynomial is a function that maps a degree configuration to an integer?? coefficient.
        For example, (2, 3) -> 4 means 4x^2y^3 (where x is the first generator and y is the second).
        """
        def __init__(self, degree_configs: list, coefficients: list):
            self.coefficient_array = []
            for n, deg_config in enumerate(degree_configs):
                self._update(deg_config, coefficients[n])

        def _update(self, deg_config: list, coefficient: int):
            target_list = self.coefficient_array
            for k, i in enumerate(deg_config):
                if i >= len(target_list):
                    target_list += [None] * (deg_config[k] - len(target_list))
                    target_list.append([])
                if k == len(deg_config) - 1:
                    target_list[i] = coefficient
                    break
                target_list = target_list[i]

        def peek(self, deg_config: list) -> int:
            target_list = self.coefficient_array
            for k, i in enumerate(deg_config):
                if k == len(deg_config) - 1:
                    return target_list[i]
                target_list = target_list[i]

        def __str__(self):
            output = ""

            return self.coefficient_array.__str__()

    class HSpace:
        """
        Space of Homogeneous Terms with the same bigrade.
        """

    def __init__(self, generators: list[str], generator_bigrades: list[tuple[int, int]]):
        self.generators = generators
        self.generator_bigrades = generator_bigrades
        self.h_spaces: list[list[BiGradedAlgebra.HSpace]] = []
        # This 2-D array is ordered in both direction to facilitate further search.
        # TODO: It might be a better practice to define a class for this

    def calculate_bigrade(self):
        pass


if __name__ == "__main__":
    BA = BiGradedAlgebra(["x", "y"], [(1, 2), (3, 4)])
    P = BA.Polynomial([[1, 4], [2, 5]], [10, 20])
    print(P)
    print(P.peek([2, 5]))