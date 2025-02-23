from numpy import ndarray, asarray


class DOArray(ndarray):
    """
    Ndarray with dictionary order.
    __eq__ is rewritten to compare the whole array.
    __hash__ calls the hash function for tuples.
    """

    def __new__(cls, input_array):
        # noinspection PyUnresolvedReferences
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
        return self.shape == other.shape and all(super().__eq__(other))

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return tuple(self).__hash__()

    # def combo(self, arr: 'DOArray') -> 'DOArray' | int:
    #     self_iter = self.__iter__()
    #     arr_iter = arr.__iter__()
    #     output = self_iter.__next__() * arr_iter.__next__()
    #     try:
    #         while True:
    #             output += self_iter.__next__() * arr_iter.__next__()
    #     except StopIteration:
    #         return output

    def __str__(self):
        return "DOArray:\n" + super().__str__()