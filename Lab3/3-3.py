import matplotlib.pyplot as plt
from LinsysSolver import *

X = [0.1, 0.5, 0.9, 1.3, 1.7, 2.1]
Y = [-2.2026, -0.19315, 0.79464, 1.5624, 2.2306, 2.8419]
N = len(X)

plt.plot(X, Y, 'ro', label='Узлы интерполяции', markersize=8)

def mnk(X, Y, num_coeffs, N):
    """
    Формирует и решает нормальную систему МНК.
    num_coeffs = степень многочлена + 1 (количество коэффициентов a0, a1, ...)
    """
    # Матрица нормальной системы: A[k][i] = Σ(x_j^(k+i))
    A = [[sum(X[j]**(k + i) for j in range(N)) for i in range(num_coeffs)] for k in range(num_coeffs)]
    x_matrix = Matrix.from_list(A)
    
    # Вектор правой части: b[k] = Σ(y_j * x_j^k)
    b = [sum(Y[j] * (X[j]**k) for j in range(N)) for k in range(num_coeffs)]
    
    a = LULinSysSolver(x_matrix, b)
    return a

def F(x, coeffs):
    """Вычисляет значение многочлена P(x) = Σ(coeffs[i] * x^i)"""
    return sum(c * x**i for i, c in enumerate(coeffs))

# Гладкая сетка для построения графиков
x_fine = [-0.1 + i/100 for i in range(250)]

a1 = mnk(X, Y, 2, N)
y1_mnk = [F(x, a1) for x in x_fine]
plt.plot(x_fine, y1_mnk, 'g-', label='МНК: степень 1', linewidth=2)

sse1 = sum((Y[i] - F(X[i], a1))**2 for i in range(N))
print(f"Коэффициенты степени 1: {a1}")
print(f"Сумма квадратов ошибок (SSE) для степени 1: {sse1:.6f}")

a2 = mnk(X, Y, 3, N)
y2_mnk = [F(x, a2) for x in x_fine]
plt.plot(x_fine, y2_mnk, 'b-', label='МНК: степень 2', linewidth=2)

sse2 = sum((Y[i] - F(X[i], a2))**2 for i in range(N))
print(f"Коэффициенты степени 2: {a2}")
print(f"Сумма квадратов ошибок (SSE) для степени 2: {sse2:.6f}")

for i in range(N):
    y1_val = F(X[i], a1)
    y2_val = F(X[i], a2)
    err1 = Y[i] - y1_val
    err2 = Y[i] - y2_val
    print(f"x={X[i]:.1f} | y={Y[i]:.4f} | P1={y1_val:.4f} | Δ1={err1:+.4f} | P2={y2_val:.4f} | Δ2={err2:+.4f}")

plt.legend()
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Приближение МНК: 1-я и 2-я степени')
plt.show()