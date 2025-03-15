"""
This files is an effort to separate linear algebra from the main program to achieve better readability.

The main setting is the following:

Consider a universe $U$. We assume:
1. Every element of $U$ is represented by a bigrade vector and a coordinate vector.
2. Elements associated with the same bigrade form a finite dimensional vector space.

We are interested in a specific kind of map on $U$. We assume:
1. For each map $f$ we have a finite set of elements on which the image of $f$ is known.
2. $f$ is a linear transformation when restricted to the vector space associated with a specific bigrade.

Additional details:
1. Every bigrades gives a "standard basis" on which every element has a representation.
2. A vector space is represented by a basis and a kernel basis. Elements spanned by the kernel basis are treated as 0.
    It's guaranteed that the basis elements are not in the kernel.
3. When comparing two elements, they are regarded as the same if the difference lies in the kernel.

Objectives:
1. Implement the structures
2. Find the matrix of a map $f$ restricted to a specific vector space.
"""

from utilities import Matrix, Vector
from sortedcontainers import SortedDict
from copy import deepcopy


class Universe:
    def __init__(self):
        self.subspaces = SortedDict([])

    def get_dimension(self, bigrade):
        return 5  # This is a test case and should be overridden.

    def get_subspace(self, bigrade) -> 'Subspace':
        if bigrade in self.subspaces:
            return self.subspaces[bigrade]
        else:
            output = self.generate_subspace(bigrade)
            self.subspaces[bigrade] = output
            return output

    def set_subspace(self, bigrade, basis, ker_basis):
        self.subspaces[bigrade] = Subspace(self, bigrade, basis, ker_basis)

    def generate_subspace(self, bigrade):
        # This is the default case. We use the standard basis and no kernel.
        return Subspace(self, bigrade, Matrix.eye(self.get_dimension(bigrade)), Matrix([]))


class HomoElement:
    """
    Class of elements

    This class will be inherited by the polynomial class, which represents elements in a bigraded algebra.

    Note: we need to take care of bigrade carefully. For a polynomial, if any part
    """
    def __init__(self, universe: Universe, bigrade: Vector):
        self.universe = universe
        self.bigrade = bigrade

    @property
    def coordinate(self) -> Vector | None:
        raise NotImplementedError

    @property
    def parent_space(self):
        assert not self.isZero()
        return self.universe.get_subspace(self.bigrade)

    def isZero(self):
        raise NotImplementedError

    def
    def __sub__(self, other: 'HomoElement') -> 'HomoElement':
        raise NotImplementedError

    def __eq__(self, other: 'Monomial'):
        return (self - other).isZero()


class Subspace:
    """
    A subspace stores a basis of itself and a kernel subspace.
    """
    def __init__(self, universe: Universe, bigrade, basis: Matrix, ker_basis: Matrix):
        """
        The basis  and ker_basis here are represented in the standard basis associated with the bigrade.
        """
        self.universe = universe
        self.bigrade = bigrade

        ker_basis_idx, basis_idx = Matrix.double_reduction(ker_basis, basis)
        self.ker_basis = Matrix([])
        for i in ker_basis_idx:
            self.ker_basis = self.ker_basis.row_join(ker_basis.col(i))
        self.basis = Matrix([])
        for i in basis_idx:
            self.basis = self.basis.row_join(basis.col(i))

    def __contains__(self, e: Monomial):
        return e.universe == self.universe and (e.isZero() or e.bigrade == self.bigrade)

    def compare(self, e1: Monomial, e2: Monomial):
        assert e1 in self and e2 in self
        return self.ker_basis.col_spans(e1.coordinate - e2.coordinate)

    def kernelContains(self, e: Monomial):
        assert e in self
        return self.ker_basis.col_spans(e.coordinate)



class LinearTransformation:
    """
    We store a linear transformation by its matrix representation. It can be defined by its image
    on a specific set of vectors.
    """

    def __init__(self, images: dict[Monomial, Monomial]):
        self.images = images

    def get_matrix(self, from_bigrade, to_bigrade):
        """
        Step 1: Take all elements with known images, generate all convex integral combination that matches the
        from_bigrade.
        Step 2: Enlarge the known elements to a basis.
        Step 3: Get the subspace corresponding to the to_bigrade and
        """