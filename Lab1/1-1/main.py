from LinsysSolver import *

a = [[21, -6, -9, -4], [-6, 20, -4, 2], [-2, -7, -20, 3], [4, 9, 6, 24]]
A = Matrix.from_list(a)

B = [127, -144, 236, 5]

print(
    f"Determinant: {Determinant(A)}\n\nSystem solution: X = {LULinSysSolver(A,B)}\n\nInverse Mtrix:\n{Inverse_matrix(A)}")
