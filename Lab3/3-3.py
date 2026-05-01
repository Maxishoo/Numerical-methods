import math
import matplotlib.pyplot as plt
from LinsysSolver import *

X = [0.1, 0.5, 0.9, 1.3, 1.7, 2.1]
Y = [-2.2026, -0.19315, 0.79464, 1.5624, 2.2306, 2.8419]
plt.plot(X, Y, 'ro', label='Узлы интерполяции', markersize=8)

N = len(X)


def mnk(X, Y, n, N):
    x = [[sum([X[j]**(k+i) for j in range(N)]) for i in range(n)] for k in range(n)]
    x_matrix = Matrix.from_list(x)# принимает список списков элементов[[a11,a12],[a21,a22]]
    b = [sum([Y[j] * (X[j]**k) for j in range(N)]) for k in range(n)]

    a = LULinSysSolver(x_matrix, b) # принимает матрицу и list
    return a

def F(x, a, n):
    ans = 0
    for i in range(n):
        ans += a[i] * x ** i
    return ans

x_fine = [-0.1 + i/100 for i in range(250)]
a2 = mnk(X,Y,3,N)
y2_mnk = [F(x,a2,3) for x in x_fine]
plt.plot(x_fine, y2_mnk, 'b-', label='приближающий многочлен степени 2', linewidth=2)
err2 = sum([(Y[i] - F(X[i], a2, 3))**2 for i in range(N)])
print(f"Коэффициенты степени 2: {a2}")
print(f"Сумма квадратов ошибок (MSE): {err2}")

a3 = mnk(X,Y,4,N)
print(a3)
y3_mnk = [F(x,a3,4) for x in x_fine]
plt.plot(x_fine, y3_mnk, 'y-', label='приближающий многочлен степени 3', linewidth=2)
err3 = sum([(Y[i] - F(X[i], a3, 4))**2 for i in range(N)])
print(f"Коэффициенты степени 3: {a3}")
print(f"Сумма квадратов ошибок (MSE): {err3}")
    


plt.legend()
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Метод наименьших квадратов')
plt.show()
