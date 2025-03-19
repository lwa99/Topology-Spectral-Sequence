# New Structure Branch Development Documentation

### Main Structure

- utilities.py

- algebra_layer.py

- spectral_sequence_layer.py

- application_layer.py

### Linear Algebra Layer

##### Universe

A universe is a set of elements. Given any bigrade, it should know the dimension of the vector space corresponding to the bigrade.

It also stores the subspaces corresponding to each bigrade.

##### Homogeneous Element

A homogeneous element is represented by a bigrade and its coordinate on the standard basis associated with the bigrade.

##### Subspace

Linear subspace of homogeneous elements with the same bigrade. Represented by a basis and a kernel basis. Homogeneous elements with coordinate spanned by the kernel basis is treated as zero.

##### Transformation

A transformation is defined by its image on some homogeneous elements. 

### Spectral Sequence Layer

##### Spectral Sequence

##### Page

Page inherits universe

##### Exponent

Exponent inherits vector. It handles the relations on generators. By modifying `__eq__`, we can treat some exponent as 1.

##### Polynomial

Polynomial is the class that handles arithmetic aspects of polynomials. Polynomials are immutable.

Every polynomial is created from monomial, which is homogeneous.

There are three kinds of polynomials. When each operation is preformed, we need to make sure which kind the output lies in:

1. Inhomogeneous Polynomials:                 bigrade is None

2. Zero Polynomial:                                        bigrade is empty

3. Homogeneous Non-zero Polynomial:    bigrade is non-trivial

Difficulties when handling polynomial is that when each polynomial is created from any operation, we need to drop off components that are actually zero. This means polynomials are classified in real time.

Information about a polynomials:

1. A sorted dictionary that map exponents to coefficients.

2. A sorted dictionary that map bigrades to coordinates. (Derived from the first one)

##### Homogeneous Polynomial

Inherits polynomial and homogeneous element

##### Module

##### Differential
