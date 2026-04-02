import math


def f(x: tuple):
    x1, x2 = x
    return 2*x1 - math.cos(x2), 2*x2 - math.exp(x1)


def jacobian_phi(x, lambd):
    x1, x2 = x
    return [
        [1 - 2*lambd, lambd * math.sin(x2)],
        [lambd * math.exp(x1), 1 - 2*lambd]
    ]


EPS = 0.000001
x_old = [0.2, 0.6]
lambd = 0.5

iters = 0
while True:
    iters += 1
    x = [0, 0]
    ff = f(x_old)
    for i in range(len(x)):
        x[i] = x_old[i]-lambd*ff[i]

    if max([abs(x_old[i]-x[i]) for i in range(len(x))]) < EPS:
        break
    x_old = x

print(f"Solution: x1 = {x_old[0]}, x2 = {x_old[1]}")
