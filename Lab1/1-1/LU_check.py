from Matrix_lib import *

A = Matrix.from_list([[10, 1, 1], [2, 10, 1], [2, 2, 10]])
L, U = LUgetter(*LUconverter(A))
print(f"RESULT:\nL:\n{L}\n\nU:\n{U}\n\nL*U:\n{L*U}\n")
