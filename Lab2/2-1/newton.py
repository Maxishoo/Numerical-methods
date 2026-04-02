import math

# x * e^x + x^2 - 1 = 0
a = 0.4
b = 0.5

def f(x):
    return x * math.exp(x) + x**2 - 1

def f_derivative(x):
    return math.exp(x) + x * math.exp(x) + 2 * x

EPS = 0.000001
MAX_ITER = 10000

x = 0.5

for i in range(MAX_ITER):
    fx = f(x)
    fpx = f_derivative(x)
    
    x_new = x - fx / fpx
    
    if abs(x_new - x) < EPS:
        print(f"Решение: {x_new}, итераций: {i+1}")
        break
    
    x = x_new
else:
    print("Превышено количество итераций")