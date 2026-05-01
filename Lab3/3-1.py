import math
import matplotlib.pyplot as plt


def f(x):
    return math.log(x) + x

def df(x, order):
    if order == 1:
        return 1 + 1/x
    elif order == 2:
        return -1/x**2
    elif order == 3:
        return 2/x**3
    elif order == 4:
        return -6/x**4
    elif order == 5:
        return 24/x**5
    else:
        return (-1)**order * math.factorial(order-1) / (x**order)

def M(a, b, order):
    return abs(df(a, order))

def omega(x, X):
    """Вычисляет ω_{n+1}(x) = ∏(x - x_i)"""
    result = 1
    for xi in X:
        result *= (x - xi)
    return result


def L(x, X, n):
    y = 0
    for i in range(n):
        tmp = f(X[i])
        for j in range(n):
            if i != j:
                tmp *= (x-X[j])/(X[i] - X[j])
        y += tmp
    return y

def split_difference_f(fi, xi, fj, xj):
    return (fi - fj) /( xi- xj)

def P(x, X, n):
    f_arr = []
    f_arr.append([f(xi) for xi in X])
    
    # Вычисляем разделённые разности высших порядков
    for order in range(1, n):
        prev = f_arr[order - 1]
        current = []
        for i in range(n - order):
            diff = (prev[i + 1] - prev[i]) / (X[i + order] - X[i])
            current.append(diff)
        f_arr.append(current)
    
    # Строим многочлен Ньютона
    result = f_arr[0][0]  # f[x₀]
    for order in range(1, n):
        term = f_arr[order][0]  # разделённая разность порядка order
        for i in range(order):  # умножаем на (x - x₀)(x - x₁)...(x - x_{order-1})
            term *= (x - X[i])
        result += term
    
    return result


def solution(X):
    X_star = 0.8
    Y = [f(i) for i in X]
    
    x_fine = [0.05 + i/100 for i in range(150)]
    y_lagrange = [L(x, X, len(X)) for x in x_fine]
    y_newton = [P(x, X, len(X)) for x in x_fine]
    y_true = [f(x) for x in x_fine]
    
    plt.plot(x_fine, y_true, 'g--', label='f(x) = ln(x) + x', alpha=0.7)
    plt.plot(X, Y, 'ro', label='Узлы интерполяции', markersize=8)

    plt.plot(x_fine, y_lagrange, 'b-', label='Многочлен Лагранжа', linewidth=2)
    plt.plot(x_fine, y_newton, '-y', label='Многочлен Ньютона', linewidth=1)
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Интерполяция многочленом Лагранжа и Ньютона')
    plt.legend()
    plt.grid(True)
    plt.show()

    eps = M(0.1, 1.3,len(X)+1)/math.factorial(len(X)+1) *abs(omega(X_star, X))
    print(f"Погрешность интерполяции по формуле: {eps}, реальная {P(X_star,X,len(X))}")


#a
solution(X=[0.1, 0.5, 0.9, 1.3])
#b
solution(X =[0.1, 0.5, 1.1, 1.3])
