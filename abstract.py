from utilities import Matrix, Vector, Prime


class AbsSpectralSequence:

    def add_relation(self, relation):
        raise NotImplementedError

    def get_basis(self, bigrade) -> tuple[Vector, ...]:
        raise NotImplementedError

    def get_dimension(self, bigrade):
        raise NotImplementedError

    def get_bigrade(self, exponent):
        raise NotImplementedError


class AbsPage:
    pass