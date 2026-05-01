import math
import matplotlib.pyplot as plt


def run_through(a: list, b: list, c: list, d: list):
    n = len(b)
    if len(a) != n or len(c) != n or len(d) != n:
        raise ValueError("Массивы разной длины")

    P = [0] * n
    Q = [0] * n

    # forward
    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]
    for i in range(1, n):
        P[i] = -c[i] / (b[i] + a[i] * P[i-1])
        Q[i] = (d[i] - a[i] * Q[i-1]) / (b[i] + a[i] * P[i-1])

    # backward
    X = [0] * n
    X[n-1] = Q[n-1]
    for i in range(n-2, -1, -1):
        X[i] = P[i] * X[i+1] + Q[i]

    return X


X_star = 0.8

X = [0.1, 0.5, 0.9, 1.3, 1.7]
Y = [-2.2026, -0.19315, 0.79464, 1.5624, 2.2306]
plt.plot(X, Y, 'ro', label='Узлы интерполяции', markersize=8)

n = len(X)
h = []
for i in range(n-1):
    h.append(X[i+1] - X[i])

# составляем систему
# размер n-2
a_system = [h[i] for i in range(1, n-1)] 
b_system = [2*(h[i-1] + h[i]) for i in range(1, n-1)]
c_system = [h[i] for i in range(1, n-1)]

d_system = [3*((Y[i+1]-Y[i])/h[i] - (Y[i] - Y[i-1])/h[i-1]) for i in range(1, n-1)]

# Добавляем нулевые граничные условия для естественного сплайна
c_full = [0] * n  # c[0] = 0, c[n-1] = 0
c_spline = run_through(a_system, b_system, c_system, d_system)
for i in range(1, n-1):
    c_full[i] = c_spline[i-1]
c = c_full

a = [Y[i] for i in range(n)]
b_coef = [0] * (n-1)
d_coef = [0] * (n-1)

for i in range(n-1):
    b_coef[i] = (Y[i+1] - Y[i])/h[i] - h[i]/3 * (2*c[i] + c[i+1])
    d_coef[i] = (c[i+1] - c[i])/(3*h[i])

def S(x_val):
    # Находим интервал, содержащий x_val
    if x_val <= X[0]:
        i = 0
    elif x_val >= X[-1]:
        i = n-2
    else:
        i = 0
        while i < n-1 and x_val > X[i+1]:
            i += 1
    
    dx = x_val - X[i]
    return a[i] + b_coef[i] * dx + c[i] * dx**2 + d_coef[i] * dx**3

x_fine = [-0.1 + i/100 for i in range(250)]
y_spline = [S(x) for x in x_fine]
plt.plot(x_fine, y_spline, 'b-', label='Сплайн', linewidth=2)
plt.legend()
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Кубический сплайн интерполяция')
plt.show()

print(f"Значение в х* = {S(X_star)}")