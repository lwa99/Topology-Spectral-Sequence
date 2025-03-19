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


class Universe:
    def __init__(self):
        self.subspaces = SortedDict([])

    def get_dimension(self, bigrade):
        # This returns the dimension of the subspace in "absolute" sense.
        return 5  # This is a test case and should be overridden.

    def get_subspace(self, bigrade) -> 'Subspace':
        if bigrade in self.subspaces:
            return self.subspaces[bigrade]
        else:
            output = self.generate_subspace(bigrade)
            self.subspaces[bigrade] = output
            return output

    def generate_subspace(self, bigrade) -> 'Subspace':
        raise NotImplementedError


class HomoElement:
    """
    Homogeneous Elements.
    """
    def __new__(cls, universe: Universe, bigrade: Vector, coordinate: Vector):
        instance = super().__new__(cls)
        instance.universe = universe
        instance.bigrade = bigrade
        instance.coordinate = coordinate
        return instance

    def isZero(self):
        """
        Since polynomials are classified in real time,
        """
        return self.bigrade == Vector([])


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

    def __contains__(self, e: HomoElement):
        return e.universe == self.universe and (e.isZero() or e.bigrade == self.bigrade)

    def kernelContains(self, vec: Vector):
        return self.ker_basis.col_spans(vec)


class LinearTransformation:
    """
    We store a linear transformation by its matrix representation. It can be defined by its image
    on a specific set of vectors.
    """

    def __init__(self, universe, images: dict[HomoElement, HomoElement]):
        self.universe: Universe = universe
        self.images = images

    def update_images(self):
        raise NotImplementedError

    def get_matrix(self, from_bigrade, to_bigrade):
        """
        Step 1: Get known elements of the same bigrade.
        Step 2: Enlarge them to a basis.
        Step 3: Join the coordinates of images.
        Step 4: Change basis back to standard
        """
        self.update_images()
        known_elements = Matrix([])
        for e in self.images.keys():
            if e.bigrade == from_bigrade:
                known_elements.row_join(e)

        # Enlarge known elements to a basis
        from_basis = self.universe.get_subspace(from_bigrade).basis
        idx1, idx2 = Matrix.double_reduction(known_elements, from_basis)

        # Now we build the temp basis and the matrix representation from temp basis to standard basis.
        temp_output = temp_basis = Matrix([])
        for i in idx1:
            temp_basis.row_join(known_elements[i])
            temp_output.row_join(self.images[known_elements[i]].coordinate)
        for i in idx2:
            temp_basis.row_join(from_basis[i])
            temp_output.row_join(Vector([0] * self.universe.get_dimension(to_bigrade)))

        return temp_output * temp_basis.inv()
