from Matrix_lib import *
import math


def SimpleItLinsysSolver(A: Matrix, B: list, Eps: float, max_iter=1000):
    n, m = A.n, A.m
    if n != m:
        raise ValueError("Число уравнений не соответствует числу неизвестных")
    if len(B) != n:
        raise ValueError("Размер вектора свободных членов не соответствует размерности матрицы")

    # 1. Проверка на диагональное преобладание
    for i in range(n):
        if A.data[i][i] == 0:
            raise ValueError(f"Нулевой элемент на диагонали в строке {i}")
        s = sum(abs(A.data[i][j]) for j in range(n) if i != j)
        if abs(A.data[i][i]) <= s:
            raise ValueError(f"Не выполнено условие диагонального преобладания в строке {i}")

    # 2. Приведение к виду x = alpha * x + beta
    beta = [B[i] / A.data[i][i] for i in range(n)]

    alpha = Matrix.zeros(n, m)
    for i in range(n):
        for j in range(m):
            if i != j:
                alpha.data[i][j] = -A.data[i][j] / A.data[i][i]

    # 3. Вычисление нормы матрицы итераций
    norm_alpha = alpha.norma()
    if norm_alpha >= 1:
        raise ValueError("Норма матрицы итераций >= 1, сходимость не гарантирована")

    k_eps = norm_alpha / (1 - norm_alpha)

    # Начальное приближение
    X = [b_val for b_val in beta]

    for iteration in range(1, max_iter + 1):
        X_new = [0.0] * n
        for i in range(n):
            sum_val = sum(alpha.data[i][j] * X[j] for j in range(m))
            X_new[i] = sum_val + beta[i]

        diff = max(abs(X[i] - X_new[i]) for i in range(n))
        if k_eps * diff < Eps:
            return X_new, iteration
        
        X = X_new[:]

    raise RuntimeError(f"Метод не сошелся за {max_iter} итераций")


def ItersEvaluation(A: Matrix, B: list, Eps: float):
    n, m = A.n, A.m
    if n != m:
        raise ValueError("Число уравнений не соответствует числу неизвестных")
    
    # Вычисление alpha и beta
    beta = [B[i] / A.data[i][i] for i in range(n)]
    alpha = Matrix.zeros(n, m)
    for i in range(n):
        for j in range(m):
            if i != j:
                alpha.data[i][j] = -A.data[i][j] / A.data[i][i]
    
    norm_alpha = alpha.norma()

    # Проверка сходимости
    if norm_alpha >= 1 or norm_alpha == 0:
        raise ValueError("Оценка невозможна: норма матрицы итераций >= 1 или равна 0")

    max_beta = max(abs(b) for b in beta)
    if max_beta == 0:
        max_beta = 1e-10
    k = (math.log10(Eps) - math.log10(max_beta) + 
            math.log10(1 - norm_alpha)) / math.log10(norm_alpha)
    return max(1, k)
