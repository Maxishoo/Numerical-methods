import math
import matplotlib.pyplot as plt

def summ(a, b):
    return [a[i]+b[i] for i in range(min(len(a), len(b)))]

def mult(a, k):
    return [k*a[i] for i in range(len(a))]

def ideal_solution(x):
    return math.sin(x) + 2 - math.sin(x) * math.log((1+math.sin(x))/(1-math.sin(x)))

def f(x, y_val, z_val):
    return [1, z_val, math.tan(x)*z_val - 2*y_val]

def runge_kutta(xa, xb, y0, z0, h):
    ans = []
    y_vec = [xa, y0, z0]
    y_vec_double = [xa, y0, z0]
    
    n_steps = int((xb - xa) / h)
    
    for i in range(n_steps + 1):
        x = y_vec[0]
        y = y_vec[1]
        
        exact_err = abs(y - ideal_solution(x))
        runge_err = abs(y - y_vec_double[1]) / 15.0
        
        ans.append((x, y, exact_err, runge_err))
        
        if i < n_steps:
            k1 = mult(f(y_vec[0], y_vec[1], y_vec[2]), h)
            k2 = mult(f(y_vec[0] + h/2, y_vec[1] + k1[1]/2, y_vec[2] + k1[2]/2), h)
            k3 = mult(f(y_vec[0] + h/2, y_vec[1] + k2[1]/2, y_vec[2] + k2[2]/2), h)
            k4 = mult(f(y_vec[0] + h, y_vec[1] + k3[1], y_vec[2] + k3[2]), h)
            sum_k = summ(k1, summ(mult(k2, 2), summ(mult(k3, 2), k4)))
            y_vec = summ(y_vec, mult(sum_k, 1/6))
            
            k1_d = mult(f(y_vec_double[0], y_vec_double[1], y_vec_double[2]), 2*h)
            k2_d = mult(f(y_vec_double[0] + h, y_vec_double[1] + k1_d[1]/2, y_vec_double[2] + k1_d[2]/2), 2*h)
            k3_d = mult(f(y_vec_double[0] + h, y_vec_double[1] + k2_d[1]/2, y_vec_double[2] + k2_d[2]/2), 2*h)
            k4_d = mult(f(y_vec_double[0] + 2*h, y_vec_double[1] + k3_d[1], y_vec_double[2] + k3_d[2]), 2*h)
            sum_k_d = summ(k1_d, summ(mult(k2_d, 2), summ(mult(k3_d, 2), k4_d)))
            y_vec_double = summ(y_vec_double, mult(sum_k_d, 1/6))
            
    return ans

def run_through(a_arr, b_arr, c_arr, d_arr):
    n = len(b_arr)
    P, Q = [0]*n, [0]*n
    P[0] = -c_arr[0] / b_arr[0]
    Q[0] = d_arr[0] / b_arr[0]
    for i in range(1, n):
        denom = b_arr[i] + a_arr[i] * P[i-1]
        P[i] = -c_arr[i] / denom
        Q[i] = (d_arr[i] - a_arr[i] * Q[i-1]) / denom
    X = [0] * n
    X[-1] = Q[-1]
    for i in range(n-2, -1, -1):
        X[i] = P[i] * X[i+1] + Q[i]
    return X

def solve_fd(a, b, h, y_a, y_b):
    n = int((b - a) / h) + 1
    x_nodes = [a + i*h for i in range(n)]
    
    def p(x): return -math.tan(x)
    def q(x): return 2
    
    A = [0] + [1/h**2 - p(a+i*h)/(2*h) for i in range(1, n-1)] + [0]
    B = [1] + [-2/h**2 + q(a+i*h) for i in range(1, n-1)] + [1]
    C = [0] + [1/h**2 + p(a+i*h)/(2*h) for i in range(1, n-1)] + [0]
    D = [y_a] + [0 for _ in range(1, n-1)] + [y_b]
    
    y_fine = run_through(A, B, C, D)
    
    h2 = 2 * h
    n2 = int((b - a) / h2) + 1
    A2 = [0] + [1/h2**2 - p(a+i*h2)/(2*h2) for i in range(1, n2-1)] + [0]
    B2 = [1] + [-2/h2**2 + q(a+i*h2) for i in range(1, n2-1)] + [1]
    C2 = [0] + [1/h2**2 + p(a+i*h2)/(2*h2) for i in range(1, n2-1)] + [0]
    D2 = [y_a] + [0 for _ in range(1, n2-1)] + [y_b]
    y_coarse = run_through(A2, B2, C2, D2)
    
    results = []
    for i in range(n):
        x = x_nodes[i]
        y = y_fine[i]
        exact_err = abs(y - ideal_solution(x))
        runge_err = abs(y - y_coarse[i//2]) / 3.0 if i % 2 == 0 else 0.0
        results.append((x, y, exact_err, runge_err))
        
    return results

# Начнем
a, b = 0, math.pi / 6
h = 0.01
y_a, y_b_target = 2, 2.5 - 0.5 * math.log(3)

z0_1, z0_2 = 1.0, 0.8
for _ in range(20):
    sol1 = runge_kutta(a, b, y_a, z0_1, h)
    sol2 = runge_kutta(a, b, y_a, z0_2, h)
    err1 = sol1[-1][1] - y_b_target
    err2 = sol2[-1][1] - y_b_target
    if abs(err2) < 1e-7: break
    z0_new = z0_2 - err2 * (z0_2 - z0_1) / (err2 - err1)
    z0_1, z0_2 = z0_2, z0_new
    
solution_shoot = runge_kutta(a, b, y_a, z0_2, h)
solution_fd = solve_fd(a, b, h, y_a, y_b_target)

print(f"{'x':<8} {'y(S)':<12} {'|y-y_точн|(S)':<14} {'|y-y_R|(S)':<12} {'y(FD)':<12} {'|y-y_точн|(FD)':<14} {'|y-y_R|(FD)':<12}")
print("-"*95)

step = max(1, len(solution_shoot) // 12)
for i in range(0, len(solution_shoot), step):
    if i < len(solution_fd):
        x_s, y_s, e_s, r_s = solution_shoot[i]
        x_f, y_f, e_f, r_f = solution_fd[i]
        y_exact = ideal_solution(x_s)
        print(f"{x_s:<8.4f} {y_s:<12.6f} {e_s:<14.2e} {r_s:<12.2e} {y_f:<12.6f} {e_f:<14.2e} {r_f:<12.2e}")
        
max_e_s = max(p[2] for p in solution_shoot)
max_r_s = max(p[3] for p in solution_shoot)
max_e_f = max(p[2] for p in solution_fd)
max_r_f = max(p[3] for p in solution_fd)

x_exact = [a + i*0.001 for i in range(int((b-a)/0.001) + 1)]
y_exact = [ideal_solution(x) for x in x_exact]

x_s = [p[0] for p in solution_shoot]
y_s = [p[1] for p in solution_shoot]
x_f = [p[0] for p in solution_fd]
y_f = [p[1] for p in solution_fd]

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(x_exact, y_exact, 'k-', linewidth=2, label='Точное решение')
plt.plot(x_f, y_f, 'bs', markersize=4, label='Конечные разности', alpha=0.7)
plt.plot(x_s, y_s, 'ro', markersize=2, label='Метод стрельбы', alpha=0.7)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Решение краевой задачи')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
err_s = [abs(y_s[i] - ideal_solution(x_s[i])) for i in range(len(x_s))]
err_f = [abs(y_f[i] - ideal_solution(x_f[i])) for i in range(len(x_f))]

plt.semilogy(x_s, err_s, 'r-', linewidth=2, label='Метод стрельбы')
plt.semilogy(x_f, err_f, 'b-', linewidth=2, label='Конечные разности')
plt.xlabel('x')
plt.ylabel('Погрешность |y_числ - y_точн|')
plt.title('Погрешности численных методов')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()