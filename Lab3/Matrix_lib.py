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

    def to_zeros(self):
        for i in range(self.n):
            for j in range(self.m):
                self.data[i][j] = 0

    def __str__(self):
        return '\n'.join(' '.join(map(str, row)) for row in self.data)


def LUconverter(A: Matrix, needcount=False):
    if A.n != A.m:
        raise ValueError("Метод применим только для квадратных матриц")

    n = A.n
    LU = A.copy()
    P = list(range(n))
    swap_count = 0

    for k in range(0, n - 1):
        max_row = k
        max_val = abs(LU.data[P[k]][k])
        for i in range(k + 1, n):
            if abs(LU.data[P[i]][k]) > max_val:
                max_val = abs(LU.data[P[i]][k])
                max_row = i
        if max_row != k:
            P[k], P[max_row] = P[max_row], P[k]
            swap_count += 1

        if LU.data[P[k]][k] == 0:
            raise ValueError("Матрица вырожденная")

        for i in range(k + 1, n):
            u = LU.data[P[i]][k] / LU.data[P[k]][k]
            LU.data[P[i]][k] = u

            for j in range(k + 1, n):
                LU.data[P[i]][j] = LU.data[P[i]][j] - u * LU.data[P[k]][j]
    if needcount:
        return P, LU, swap_count
    else:
        return P, LU


def LUgetter(P: list, LU: Matrix):
    U = Matrix.zeros(LU.n, LU.m)
    L = Matrix.zeros(LU.n, LU.m)
    for i in range(LU.n):
        for j in range(LU.m):
            if i <= j:
                if i == j:
                    L.data[i][j] = 1
                U.data[i][j] = LU.data[P[i]][j]
            else:
                L.data[i][j] = LU.data[P[i]][j]
    return L, U
