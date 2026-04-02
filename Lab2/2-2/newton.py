import math

class Matrix:
    def __init__(self, data, n, m):
        self.data = data
        self.n = n
        self.m = m

    @classmethod
    def from_list(cls, matrix):
        n = len(matrix)
        m = len(matrix[0]) if n > 0 else 0
        return cls([row[:] for row in matrix], n, m)

    @classmethod
    def zeros(cls, n, m):
        return cls([[0 for _ in range(m)] for _ in range(n)], n, m)

    def copy(self):
        new_matrix = [row[:] for row in self.data]
        return Matrix(new_matrix, self.n, self.m)

    def __getitem__(self, idx):
        return self.data[idx]

    def __setitem__(self, idx, value):
        self.data[idx] = value

    def __mul__(self, other):
        if self.m != other.n:
            raise ValueError("Несовместимые размеры матриц")
        ans = []
        for i in range(self.n):
            cur = []
            for j in range(other.m):
                sum_val = 0
                for k in range(self.m):
                    sum_val += self.data[i][k] * other.data[k][j]
                cur.append(sum_val)
            ans.append(cur)
        return Matrix(ans, self.n, other.m)
    
    def __add__(self, other):
        if self.n != other.n or self.m != other.m:
            raise ValueError("Несовместимые размеры матриц")
        ans = []
        for i in range(self.n):
            tmp = []
            for j in range(self.m):
                tmp.append(self.data[i][j] + other.data[i][j])
            ans.append(tmp)
        return Matrix(ans, self.n, self.m)

    def norma(self):
        ans = 0
        for i in range(self.n):
            s = 0
            for j in range(self.m):
                s += abs(self.data[i][j])
            ans = max(ans, s)
        return ans

    def __str__(self):
        return '\n'.join(' '.join(map(str, row)) for row in self.data)


def SeidelLinsysSolver(A: Matrix, B: list, Eps: float, max_iter=1000):
    n, m = A.n, A.m
    
    X = [0.0] * n

    for iteration in range(1, max_iter + 1):
        max_diff = 0.0
        for i in range(n):
            sum_val = sum(A.data[i][j] * X[j] for j in range(n) if j != i)
            x_new = (B[i] - sum_val) / A.data[i][i]
            
            diff = abs(x_new - X[i])
            if diff > max_diff:
                max_diff = diff
            
            X[i] = x_new

        if max_diff < Eps:
            return X, iteration

    raise RuntimeError(f"Метод не сошелся за {max_iter} итераций")


def f(x):
    x1, x2 = x
    return [2*x1 - math.cos(x2), 2*x2 - math.exp(x1)]


def jacobian_f(x):
    x1, x2 = x
    return [
        [2, math.sin(x2)],
        [-math.exp(x1), 2]
    ]

EPS = 0.000001
MAX_ITER = 1000

x = [0.5, 0.5]

for i in range(MAX_ITER):
    fx = f(x)
    
    residual = max(abs(val) for val in fx)
    
    J = Matrix.from_list(jacobian_f(x))
    
    # Решаем систему J · Δx = -f(x)
    neg_fx = [-val for val in fx]
    dx, iters_inner = SeidelLinsysSolver(J, neg_fx, EPS/10)

    x_new = [x[i] + dx[i] for i in range(len(x))]
    
    if max(abs(x_new[i] - x[i]) for i in range(len(x))) < EPS:
        print(f"Решение: x1 = {x_new[0]}, x2 = {x_new[1]} Итераций Ньютона: {i+1}")
        break
    
    x = x_new
else:
    print("Превышено количество итераций")