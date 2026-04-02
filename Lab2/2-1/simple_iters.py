import math

# x * e^x + x^2 - 1 = 0
a = 0.4
b = 0.5

EPS = 0.000001


def f(x):
    return x * math.exp(x) + x**2 - 1


def f_derivative(x):
    return math.exp(x) + x * math.exp(x) + 2 * x


x_old = 0.4

lambd = 0.5
print(f"try lambda = {abs(1 - lambd * f_derivative(x_old))}")

iters = 0
while True:
    iters += 1
    x = x_old - lambd * f(x_old)
    if abs(x_old - x) < EPS:
        break
    x_old = x

print(f"Solution: {x}, iters: {iters}")
