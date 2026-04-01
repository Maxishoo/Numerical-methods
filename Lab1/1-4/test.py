from RotateLib import *

def check(A_original: Matrix, eigenvalues, eigenvectors, eps=1e-6):
    n = A_original.n
    results = []
    max_error = 0
    
    for i in range(n):
        λ = eigenvalues[i]
        v = eigenvectors[i]
        
        # 1. Вычисляем A · v
        Av = [0] * n
        for row in range(n):
            for col in range(n):
                Av[row] += A_original.data[row][col] * v[col]
        
        # 2. Вычисляем λ · v
        λv = [λ * v[j] for j in range(n)]
        
        # 3. ||Av - λv||
        error = math.sqrt(sum((Av[j] - λv[j])**2 for j in range(n)))
        max_error = max(max_error, error)
        
        results.append(error < eps)
        print(f"Вектор #{i+1}: λ={λ:.6f}, {'yes' if error < eps else 'no'}")
    print(f"\nМаксимальная ошибка: {max_error:.2e}")
    return all(results), max_error


a = [[8,-3,9],[-3,8,-2],[9,-2,-8]]
A = Matrix.from_list(a)
A_original = A.copy()

eigenvalues, eigenvectors, iters = RotateMethod(A, Eps=1e-8)
print("=== ПРОВЕРКА: A·v = λ·v ===")
all_ok, max_err = check(A_original, eigenvalues, eigenvectors)