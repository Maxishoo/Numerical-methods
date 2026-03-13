from Matrix_lib import *


def SeidelLinsysSolver(A: Matrix, B: list, Eps: float, max_iter=1000):
    n, m = A.n, A.m
    if n != m:
        raise ValueError("Число уравнений не соответствует числу неизвестных")
    if len(B) != n:
        raise ValueError("Размер вектора свободных членов не соответствует размерности матрицы")

    # Проверка диагонального преобладания
    for i in range(n):
        if A.data[i][i] == 0:
            raise ValueError(f"Нулевой элемент на диагонали в строке {i}")
        s = sum(abs(A.data[i][j]) for j in range(n) if i != j)
        if abs(A.data[i][i]) <= s:
            raise ValueError(f"Нет диагонального преобладания в строке {i}")

    X = [0.0] * n

    for iteration in range(1, max_iter + 1):
        max_diff = 0.0
        for i in range(n):
            sum_val = sum(A.data[i][j] * X[j] for j in range(n) if j != i)
            x_new = (B[i] - sum_val) / A.data[i][i]
            
            diff = abs(x_new - X[i])
            if diff > max_diff:
                max_diff = diff
            
            X[i] = x_new

        if max_diff < Eps:
            return X, iteration

    raise RuntimeError(f"Метод не сошелся за {max_iter} итераций")
