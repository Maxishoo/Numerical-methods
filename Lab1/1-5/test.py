from QRLib import *

a = [[1,2,5],[-8,0,-6],[7,-9,-7]]
A = Matrix.from_list(a)

print(QREigenvalues(A, 0.00001))
# https://mathforyou.net/online/matrices/eigen/