from Matrix_lib import *
import math

def RotateMethod(A: Matrix, Eps: float):
    if not is_symmetric(A):
        raise ValueError("Матрица должна быть симметричной")
    
    n = A.n
    Vectors = Matrix.singular(n, n)
    iters = 0
    
    while True:
        # Поиск максимального по модулю внедиагонального элемента
        mx = 0
        l, p = 0, 1
        for i in range(n):
            for j in range(i):
                if abs(A.data[i][j]) > mx:
                    mx = abs(A.data[i][j])
                    l, p = i, j
        
        # Критерий остановки
        if mx < Eps:
            break
        
        # Вычисление угла вращения
        diag_diff = A.data[l][l] - A.data[p][p]
        if abs(diag_diff) < 1e-12:
            fi = math.pi / 4 if A.data[l][p] > 0 else -math.pi / 4
        else:
            fi = 0.5 * math.atan2(2 * A.data[l][p], diag_diff)
        
        # Матрица вращения
        U = Matrix.singular(n, n)
        c, s = math.cos(fi), math.sin(fi)
        U.data[l][l] = c
        U.data[l][p] = -s
        U.data[p][l] = s
        U.data[p][p] = c
        
        # Применение вращения
        A = U.transpose() * A * U
        Vectors = Vectors * U
        iters += 1
    
    eigenvalues = [A.data[i][i] for i in range(n)]
    eigenvectors = Vectors.transpose().data
    return eigenvalues, eigenvectors, iters

def is_symmetric(A: Matrix, eps=1e-10):
    if A.n != A.m:
        return False
    for i in range(A.n):
        for j in range(i+1, A.m):
            if abs(A.data[i][j] - A.data[j][i]) > eps:
                return False
    return True