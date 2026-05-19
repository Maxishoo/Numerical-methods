import math
import matplotlib.pyplot as plt

# --- Вспомогательные функции ---
def summ(a, b):
    return [a[i]+b[i] for i in range(min(len(a), len(b)))]

def mult(a, k):
    return [k*a[i] for i in range(len(a))]

def ideal_solution(x):
    # Точное решение
    sin_x = math.sin(x)
    # Защита от логарифма отрицательного числа или нуля
    arg = (1 + sin_x) / (1 - sin_x)
    if arg <= 0: return 0 
    return sin_x + 2 - sin_x * math.log(arg)

def f_ode(x, y_val, z_val):
    # Система: y' = z, z' = tan(x)*z - 2*y
    return [1, z_val, math.tan(x)*z_val - 2*y_val]

def runge_kutta_step(x, y, z, h):
    """Один шаг классического RK4"""
    k1 = f_ode(x, y, z)
    k2 = f_ode(x + h/2, y + h*k1[1]/2, z + h*k1[2]/2)
    k3 = f_ode(x + h/2, y + h*k2[1]/2, z + h*k2[2]/2)
    k4 = f_ode(x + h, y + h*k3[1], z + h*k3[2])
    
    dy = h/6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    dz = h/6 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    return x + h, y + dy, z + dz

def runge_kutta_full(xa, xb, y0, z0, h):
    """Полное решение задачи Коши с оценкой Рунге"""
    ans = []
    x, y, z = xa, y0, z0
    n_steps = int((xb - xa) / h)
    
    # Для оценки Рунге
    y_prev, z_prev = y0, z0
    x_prev = xa
    
    for i in range(n_steps + 1):
        exact_err = abs(y - ideal_solution(x))
        
        # Оценка Рунге (правильная для метода 4-го порядка)
        if i > 0:
            _, y_2h, _ = runge_kutta_step(x_prev, y_prev, z_prev, 2*h)
            runge_err = abs(y - y_2h) / 15.0  # 2^4 - 1 = 15
        else:
            runge_err = 0.0
            
        ans.append((x, y, exact_err, runge_err))
        
        if i < n_steps:
            x_prev, y_prev, z_prev = x, y, z
            x, y, z = runge_kutta_step(x, y, z, h)
            
    return ans

def runge_kutta_simple(xa, xb, y0, z0, h):
    """Простое решение без сохранения промежуточных точек (для метода стрельбы)"""
    x, y, z = xa, y0, z0
    n_steps = int((xb - xa) / h)
    
    for _ in range(n_steps):
        x, y, z = runge_kutta_step(x, y, z, h)
    
    return y  # возвращаем только y(b)

def run_through(a_arr, b_arr, c_arr, d_arr):
    """Простая прогонка для трехдиагональной матрицы"""
    n = len(b_arr)
    if n == 0: return []
    P, Q = [0]*n, [0]*n
    
    # Прямой ход
    for i in range(n):
        if i == 0:
            denom = b_arr[i]
            if abs(denom) < 1e-12: continue
            P[i] = -c_arr[i] / denom
            Q[i] = d_arr[i] / denom
        else:
            denom = b_arr[i] + a_arr[i] * P[i-1]
            if abs(denom) < 1e-12: continue
            P[i] = -c_arr[i] / denom
            Q[i] = (d_arr[i] - a_arr[i] * Q[i-1]) / denom
            
    # Обратный ход
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
    
    # Коэффициенты разностной схемы
    A = [0]*n; B = [0]*n; C = [0]*n; D = [0]*n
    
    # Граничные условия
    B[0] = 1; D[0] = y_a
    B[-1] = 1; D[-1] = y_b
    
    # Внутренние узлы
    for i in range(1, n-1):
        xi = x_nodes[i]
        A[i] = 1/h**2 - p(xi)/(2*h)
        B[i] = -2/h**2 + q(xi)
        C[i] = 1/h**2 + p(xi)/(2*h)
        D[i] = 0
        
    y_fine = run_through(A, B, C, D)
    
    # Оценка Рунге для КР
    h2 = 2 * h
    n2 = int((b - a) / h2) + 1
    x_coarse = [a + i*h2 for i in range(n2)]
    
    A2 = [0]*n2; B2 = [0]*n2; C2 = [0]*n2; D2 = [0]*n2
    B2[0] = 1; D2[0] = y_a
    B2[-1] = 1; D2[-1] = y_b
    
    for i in range(1, n2-1):
        xi = x_coarse[i]
        A2[i] = 1/h2**2 - p(xi)/(2*h2)
        B2[i] = -2/h2**2 + q(xi)
        C2[i] = 1/h2**2 + p(xi)/(2*h2)
        
    y_coarse = run_through(A2, B2, C2, D2)
    
    results = []
    for i in range(n):
        x = x_nodes[i]
        y = y_fine[i]
        exact_err = abs(y - ideal_solution(x))
        
        runge_err = 0.0
        if i % 2 == 0 and i//2 < len(y_coarse):
             runge_err = abs(y - y_coarse[i//2]) / 3.0  # 2^2 - 1 = 3 для метода 2-го порядка
             
        results.append((x, y, exact_err, runge_err))
        
    return results

# ================= MAIN =================
a, b = 0, math.pi / 6
h = 0.01
y_a, y_b_target = 2, 2.5 - 0.5 * math.log(3)

print(f"Точное значение y(b) = {y_b_target:.10f}")
print(f"Точное решение в точке a: y(a)={ideal_solution(a):.10f}")
print(f"Точное решение в точке b: y(b)={ideal_solution(b):.10f}\n")

# --- Метод стрельбы с улучшенным алгоритмом ---
# Используем метод Ньютона с аналитическим решением вариационного уравнения
# Но для простоты используем улучшенный метод секущих с демпфированием

# Оцениваем начальное приближение из граничных условий
z_guess = (y_b_target - y_a) / (b - a)  # линейное приближение

z_prev = z_guess
y_prev = runge_kutta_simple(a, b, y_a, z_prev, h)
err_prev = y_prev - y_b_target

print(f"Начальное приближение z0 = {z_prev:.6f}, ошибка в b: {err_prev:.6e}")

# Пробуем второе приближение
z_curr = z_prev * 1.1  # небольшое отклонение
y_curr = runge_kutta_simple(a, b, y_a, z_curr, h)
err_curr = y_curr - y_b_target

print(f"Второе приближение z0 = {z_curr:.6f}, ошибка в b: {err_curr:.6e}")

# Метод секущих с демпфированием
max_iter = 50
tol = 1e-10
damping = 1.0  # коэффициент демпфирования

for iteration in range(max_iter):
    if abs(err_curr) < tol:
        print(f"\nСходимость достигнута за {iteration+1} итераций")
        break
    
    # Вычисляем новое приближение
    derivative = (err_curr - err_prev) / (z_curr - z_prev)
    if abs(derivative) < 1e-15:
        print("Производная близка к нулю, используем другое приближение")
        z_new = z_curr * 1.5
    else:
        # Метод Ньютона с демпфированием
        z_new = z_curr - damping * err_curr / derivative
    
    # Проверяем на выход за допустимые пределы
    if abs(z_new) > 1e6:
        print("Решение расходится, изменяем начальное приближение")
        z_new = z_guess * (1 + 0.1 * iteration)
    
    # Вычисляем новую ошибку
    y_new = runge_kutta_simple(a, b, y_a, z_new, h)
    err_new = y_new - y_b_target
    
    print(f"Итерация {iteration+1}: z0 = {z_new:.10f}, ошибка = {err_new:.6e}")
    
    # Проверка на расходимость
    if abs(err_new) > abs(err_curr) * 2:
        damping *= 0.5  # уменьшаем шаг
        print(f"  Демпфирование уменьшено до {damping:.3f}")
        continue
    
    # Обновляем значения для следующей итерации
    z_prev, err_prev = z_curr, err_curr
    z_curr, err_curr = z_new, err_new
    damping = min(1.0, damping * 1.2)  # постепенно увеличиваем демпфирование

# Финальное решение с найденным z0
z_final = z_curr
print(f"\nФинальное значение z(a) = {z_final:.10f}")

# Получаем полное решение для построения графиков
solution_shoot = runge_kutta_full(a, b, y_a, z_final, h)
solution_fd = solve_fd(a, b, h, y_a, y_b_target)

# --- Вывод таблицы результатов ---
print(f"\n{'x':<8} {'y(S)':<10} {'Err(S)':<10} {'R(S)':<10} || {'y(FD)':<10} {'Err(FD)':<10} {'R(FD)':<10}")
print("-" * 85)

step = max(1, len(solution_shoot) // 10)
for i in range(0, len(solution_shoot), step):
    if i < len(solution_fd):
        x_s, y_s, e_s, r_s = solution_shoot[i]
        x_f, y_f, e_f, r_f = solution_fd[i]
        print(f"{x_s:<8.4f} {y_s:<10.6f} {e_s:<10.2e} {r_s:<10.2e} || {y_f:<10.6f} {e_f:<10.2e} {r_f:<10.2e}")

# --- Построение графиков ---
x_exact = [a + i*0.001 for i in range(int((b-a)/0.001) + 1)]
y_exact = [ideal_solution(x) for x in x_exact]

x_s = [p[0] for p in solution_shoot]
y_s = [p[1] for p in solution_shoot]
x_f = [p[0] for p in solution_fd]
y_f = [p[1] for p in solution_fd]

plt.figure(figsize=(14, 6))

# График решений
plt.subplot(1, 2, 1)
plt.plot(x_exact, y_exact, 'k-', linewidth=2, label='Точное решение')
plt.plot(x_f, y_f, 'bs', markersize=4, label=f'Конечные разности (h={h})', alpha=0.7)
plt.plot(x_s, y_s, 'ro', markersize=3, label='Метод стрельбы (RK4)', alpha=0.7)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Сравнение методов')
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3)

# График погрешностей
plt.subplot(1, 2, 2)
err_s = [p[2] for p in solution_shoot]
err_f = [p[2] for p in solution_fd]

plt.semilogy(x_s, err_s, 'r-', linewidth=1.5, label='Погрешность стрельбы (RK4)')
plt.semilogy(x_f, err_f, 'b-', linewidth=1.5, label='Погрешность сеточного метода')
plt.xlabel('x')
plt.ylabel('|y_числ - y_точн|')
plt.title('Абсолютная погрешность (Log Scale)')
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3, which="both")

plt.tight_layout()
plt.show()