import math
from Matrix_lib import *

def sign(n):
    if n > 0:
        return 1
    if n < 0:
        return -1
    return 0

def QREigenvalues(A: Matrix, Eps: float, max_iters=1000):
    n = A.n
    if A.n != A.m:
        raise ValueError("Матрица должна быть квадратной")
    
    A = A.copy()
    iters = 0
    
    while iters < max_iters:
        Q_total = Matrix.singular(n, n)
        
        for j in range(n - 1):
            norm = math.sqrt(sum(A.data[i][j]**2 for i in range(j, n)))
            
            if norm < Eps:
                continue
            
            # Вектор Хаусхолдера
            v = [0] * n
            s = sign(A.data[j][j])
            v[j] = A.data[j][j] + s * norm
            
            for i in range(j + 1, n):
                v[i] = A.data[i][j]
            
            # Норма вектора v
            v_norm_sq = sum(x**2 for x in v)
            if v_norm_sq < Eps:
                continue
            
            # Матрица Хаусхолдера H = E - 2*v*v^T / (v^T*v)
            H = Matrix.singular(n, n)
            for i in range(n):
                for k in range(n):
                    H.data[i][k] -= 2 * v[i] * v[k] / v_norm_sq
            
            # Применяем отражение: A = H * A, Q = Q * H
            A = H * A
            Q_total = Q_total * H
        
        R = A
        
        # === Шаг 2: A = R * Q ===
        A = R * Q_total
        iters += 1
        
        # === Шаг 3: Проверка сходимости (внедиагональные элементы) ===
        max_offdiag = 0
        for i in range(n):
            for j in range(i):  # Нижняя треугольная часть
                max_offdiag = max(max_offdiag, abs(A.data[i][j]))
        
        if max_offdiag < Eps:
            break
    
    # === Шаг 4: Извлечение собственных значений ===
    eigenvalues = []
    i = 0
    while i < n:
        if i == n - 1:
            # Последний элемент — вещественное собственное значение
            eigenvalues.append(A.data[i][i])
            i += 1
        elif abs(A.data[i+1][i]) < Eps:
            # Диагональный элемент — вещественное собственное значение
            eigenvalues.append(A.data[i][i])
            i += 1
        else:
            # 2×2 блок — комплексное собственное значение
            a = A.data[i][i]
            b = A.data[i][i+1]
            c = A.data[i+1][i]
            d = A.data[i+1][i+1]
            
            # Характеристическое уравнение: λ² - (a+d)λ + (ad-bc) = 0
            trace = a + d
            det = a * d - b * c
            discriminant = trace**2 - 4 * det
            
            if discriminant >= 0:
                # Вещественные корни
                sqrt_d = math.sqrt(discriminant)
                eigenvalues.append((trace + sqrt_d) / 2)
                eigenvalues.append((trace - sqrt_d) / 2)
            else:
                # Комплексные корни
                real_part = trace / 2
                imag_part = math.sqrt(-discriminant) / 2
                eigenvalues.append((real_part, imag_part))
                eigenvalues.append((real_part, -imag_part))
            
            i += 2
    
    return eigenvalues, iters