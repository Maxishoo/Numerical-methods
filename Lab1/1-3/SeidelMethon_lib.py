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
    
    norm_C = 0.0
    norm_alpha = 0.0
    
    for i in range(n):
        row_C = 0.0
        row_alpha = 0.0
        
        for j in range(n):
            if i != j:
                alpha_ij = abs(A.data[i][j] / A.data[i][i])
                row_alpha += alpha_ij
                
                if j > i:
                    row_C += alpha_ij
        
        norm_C = max(norm_C, row_C)
        norm_alpha = max(norm_alpha, row_alpha)
    
    if norm_alpha >= 1:
        raise ValueError(f"Норма матрицы итераций = {norm_alpha} >= 1, нет гарантированной сходимости")
    
    k_eps = norm_C / (1 - norm_alpha)

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

        if max_diff*k_eps < Eps:
            return X, iteration

    raise RuntimeError(f"Метод не сошелся за {max_iter} итераций")

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
    norm_C = 0.0
    norm_alpha = 0.0
    
    for i in range(n):
        row_C = 0.0
        row_alpha = 0.0
        
        for j in range(n):
            if i != j:
                alpha_ij = abs(A.data[i][j] / A.data[i][i])
                row_alpha += alpha_ij
                
                if j > i:
                    row_C += alpha_ij
        
        norm_C = max(norm_C, row_C)
        norm_alpha = max(norm_alpha, row_alpha)
    
    if norm_alpha >= 1:
        raise ValueError(f"Норма матрицы итераций = {norm_alpha} >= 1")
    
    k_eps = norm_C / (1 - norm_alpha)
    X = [0.0] * n
    X_prev = [0.0] * n
    for iteration in range(1, max_iter + 1):
        X_prev = X[:]
        max_diff = 0.0

        for i in range(n):
            sum_val = sum(A.data[i][j] * X[j] for j in range(n) if j != i)
            x_new = (B[i] - sum_val) / A.data[i][i]
            max_diff = max(max_diff, abs(x_new - X[i]))
            X[i] = x_new
        error_estimate = k_eps * max_diff
        
        if error_estimate < Eps:
            print(f"Итерация {iteration}: оценка погрешности = {error_estimate:.2e}")
            return X, iteration

    raise RuntimeError(f"Метод не сошелся за {max_iter} итераций")