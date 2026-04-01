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

    @classmethod
    def singular(cls, n, m):
        A = [[0]*m for _ in range(n)]
        for i in range(min(n, m)):
            A[i][i] = 1
        return cls(A, n, m)

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

    def transpose(self):
        ans = Matrix.zeros(self.m, self.n)
        for i in range(self.n):
            for j in range(self.m):
                ans.data[j][i] = self.data[i][j]
        return ans

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


def supportLinSysSolver(P: list, LU: Matrix, B: list):
    n = LU.n

    PB = [B[P[i]] for i in range(n)]
    z = []
    for i in range(0, n):
        z.append(PB[i])
        for j in range(0, i):
            z[i] -= LU.data[P[i]][j] * z[j]

    # PUx = Pz
    ans = [0] * n
    for i in range(n-1, -1, -1):
        tmp = z[i]
        for j in range(i+1, n):
            tmp -= LU.data[P[i]][j] * ans[j]
        tmp /= LU.data[P[i]][i]
        ans[i] = tmp

    return ans


def LULinSysSolver(A: Matrix, B: list):
    P, LU = LUconverter(A)
    return supportLinSysSolver(P, LU, B)


def Inverse_matrix(A: Matrix):
    P, LU = LUconverter(A)
    n = LU.n

    e = [0] * n
    ansSt = []
    for i in range(n):
        e[i] = 1
        ansSt.append(supportLinSysSolver(P, LU, e))
        e[i] = 0

    res = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            res.data[i][j] = ansSt[j][i]
    return res


def Determinant(A: Matrix):
    if A.n != A.m:
        raise ValueError("Матрица не является квадратной")

    P, LU, swap_count = LUconverter(A, True)
    det = 1.0
    for i in range(A.n):
        det *= LU.data[P[i]][i]
    return det * ((-1)**swap_count)
