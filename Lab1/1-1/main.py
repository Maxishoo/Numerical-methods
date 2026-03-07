from LinsysSolver import *

a = [[-5, -1, -3, -1], [-2, 0, 8, -4], [-7, -2, 2, -2], [2, -4, -4, 4]]
A = Matrix.from_list(a)

B = [18, -12, 6, -12]

print(
    f"Determinant: {Determinant(A)}\n\nSystem solution: X = {LinSysSolver(A,B)}")
