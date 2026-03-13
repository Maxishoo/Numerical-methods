from SeidelMethon_lib import *
from SimpleIterations_lib import*


a = [[21, -6, -9, -4], [-6, 20, -4, 2], [-2, -7, -20, 3], [4, 9, 6, 24]]
A = Matrix.from_list(a)

B = [127, -144, 236, 5]

EPS = 0.0001

print(f"Iters need: {ItersEvaluation(A,B,EPS)}\n\nSimpleIters: {SimpleItLinsysSolver(A,B,EPS)}\n\nSeidel: {SeidelLinsysSolver(A,B,EPS)}\n\n")