
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

# solution
a = [0, 2, -9, -4, 7]
b = [18, -9, 21, -10, 12]
c = [-9, -4, -8, 5, 0]
d = [-81, 71, -39, 64, 3]

print(f"Solution: {run_through(a,b,c,d)}")
