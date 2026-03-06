from LinsysSolver import *

a = [[-5, -1, -3,-1],[-2,0,8,-4],[-7,-2,2,-2],[2,-4,-4,4]]
A = Matrix.from_list(a)

B = [18,-12,6,-12]

print(f"Determinant: {Determinant(A)}\n\nSystem solution: X = {LinSysSolver(A,B)}")


a = [[-2, -1, -9,-5],[-4,4,-2,6],[0,5,7,-4],[0,9,7,7]]
A = Matrix.from_list(a)

B = [93,16,-80,-119]

print(f"Determinant: {Determinant(A)}\n\nSystem solution: X = {LinSysSolver(A,B)}")