from Matrix_lib import *

def supportLinSysSolver(P: list, LU: Matrix, B: list):
    n = LU.n
    
    PB = [B[P[i]] for i in range(n)]
    z = []
    for i in range(0, n):
        z.append(PB[i])
        for j in range(0, i):
            z[i] -= LU.data[P[i]][j] * z[j]
    
    # PUx = Pz
    ans = [0] * n  
    for i in range(n-1, -1, -1):
        tmp = z[i]
        for j in range(i+1, n):
            tmp -= LU.data[P[i]][j] * ans[j]
        tmp /= LU.data[P[i]][i]
        ans[i] = tmp
    
    return ans

def LinSysSolver(A: Matrix, B: list):
    P, LU = LUconverter(A)
    return supportLinSysSolver(P, LU, B)

def Inverse_matrix(A: Matrix):
    P, LU = LUconverter(A)
    n = LU.n
    
    e = [0] * n
    ansSt = []
    for i in range(n):
        e[i] = 1
        ansSt.append(supportLinSysSolver(P, LU, e))
        e[i] = 0
    
    res = Matrix.zeros(n,n)
    for i in range(n):
        for j in range(n):
            res.data[i][j] = ansSt[j][i]
    return res

def Determinant(A: Matrix):
    if A.n != A.m:
        raise ValueError("Матрица не является квадратной")
    
    P, LU, swap_count = LUconverter(A, True)
    det = 1.0
    for i in range(A.n):
        det *= LU.data[P[i]][i]
    return det * ((-1)**swap_count)