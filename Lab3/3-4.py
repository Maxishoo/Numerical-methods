import math
import matplotlib.pyplot as plt

X = [0.0, 1.0, 2.0, 3.0, 4.0]
Y = [0.0, 2.0, 3.4142, 4.7321, 6.0]
plt.plot(X, Y, 'ro', label='Узлы', markersize=8)

X_star = 2.0

def df(x):
    i = 0
    while i < len(X) - 1 and x > X[i]:
        i += 1
    
    if i == 0:
        i = 1
    if i >= len(X) - 1:
        i = len(X) - 2
    
    x_i = X[i-1]
    x_i1 = X[i]
    x_i2 = X[i+1]
    y_i = Y[i-1]
    y_i1 = Y[i]
    y_i2 = Y[i+1]

    term1 = (y_i1 - y_i) / (x_i1 - x_i)
    term2_num = (y_i2 - y_i1) / (x_i2 - x_i1) - (y_i1 - y_i) / (x_i1 - x_i)
    term2 = term2_num / (x_i2 - x_i) * (2*x - x_i - x_i1)
    
    return term1 + term2

def d2f(x):
    i = 0
    while i < len(X) - 1 and x > X[i]:
        i += 1
    
    if i == 0:
        i = 1
    if i >= len(X) - 1:
        i = len(X) - 2
    
    x_i = X[i-1]
    x_i1 = X[i]
    x_i2 = X[i+1]
    y_i = Y[i-1]
    y_i1 = Y[i]
    y_i2 = Y[i+1]
    
    numerator = (y_i2 - y_i1) / (x_i2 - x_i1) - (y_i1 - y_i) / (x_i1 - x_i)
    denominator = x_i2 - x_i
    
    return 2 * numerator / denominator

first_derivative = df(X_star)
second_derivative = d2f(X_star)

print(f"Первая производная в точке x = {X_star}: {first_derivative}")
print(f"Вторая производная в точке x = {X_star}: {second_derivative}")

plt.axvline(x=X_star, color='g', linestyle='--', label=f'x* = {X_star}')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Численное дифференцирование')
plt.legend()
plt.grid(True)
plt.show()