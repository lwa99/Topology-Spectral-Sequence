from utilities import Matrix, pprint

a = Matrix([[1, 2, 3]])
b = Matrix([4])
pprint(a.row_join(b))
pprint(a)