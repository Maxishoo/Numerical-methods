import math
import matplotlib.pyplot as plt

#прогонка
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


def summ(a, b):
    return [a[i]+b[i] for i in range(min(len(a), len(b)))]


def mult(a, k):
    return [k*a[i] for i in range(len(a))]


# x in [2,3]
h = 0.1
# (x**2-1)y''-2xy'+2y=0    y(2) = 7  y`(2) = 5
# y'' = (2xy' - 2y)/(x^2-1)
# Система: x'=1, y'=z, z'=(2xz-2y)/(x^2-1)


def ideal_ans(x):
    return x**2 + x + 1

def f(x, y_val, z_val):
    return [
        1,
        z_val,
        (2*x*z_val - 2*y_val)/(x**2 - 1)
    ]


# Метод эйлера 1го порядка точности
def eiler(x0, y0, z0, h):
    ans = [(x0, y0)]
    y_vec = [x0, y0, z0]
    y_vec_double = [x0, y0, z0]
    
    for i in range(10):
        dy = f(y_vec[0], y_vec[1], y_vec[2])
        y_vec = summ(y_vec, mult(dy, h))
        ans.append((y_vec[0], y_vec[1]))
        
        dy_double = f(y_vec_double[0], y_vec_double[1], y_vec_double[2])
        y_vec_double = summ(y_vec_double, mult(dy_double, 2*h))
        
        runge_romberg_err = abs(y_vec[1] - y_vec_double[1])
        exact_err = abs(y_vec[1] - ideal_ans(y_vec[0]))
        
        print(f"In x = {y_vec[0]:.2f} | exact_error = {exact_err:.6f} | Runge-Romberg = {runge_romberg_err:.6f}")
    
    return ans


print("=== Euler Method ===")
eiler_result = eiler(2, 7, 5, h)


# Метод Рунге-Кутты 4го порядка точности
def runge_kutta(x0, y0, z0, h):
    ans = [(x0, y0)]
    y_vec = [x0, y0, z0]
    y_vec_double = [x0, y0, z0]
    
    for i in range(10):
        k1 = mult(f(y_vec[0], y_vec[1], y_vec[2]), h)
        k2 = mult(f(y_vec[0] + h/2, y_vec[1] + k1[1]/2, y_vec[2] + k1[2]/2), h)
        k3 = mult(f(y_vec[0] + h/2, y_vec[1] + k2[1]/2, y_vec[2] + k2[2]/2), h)
        k4 = mult(f(y_vec[0] + h, y_vec[1] + k3[1], y_vec[2] + k3[2]), h)

        sum_k = summ(k1, summ(mult(k2, 2), summ(mult(k3, 2), k4)))
        y_vec = summ(y_vec, mult(sum_k, 1/6))
        ans.append((y_vec[0], y_vec[1]))
        
        k1_d = mult(f(y_vec_double[0], y_vec_double[1], y_vec_double[2]), 2*h)
        k2_d = mult(f(y_vec_double[0] + h, y_vec_double[1] + k1_d[1]/2, y_vec_double[2] + k1_d[2]/2), 2*h)
        k3_d = mult(f(y_vec_double[0] + h, y_vec_double[1] + k2_d[1]/2, y_vec_double[2] + k2_d[2]/2), 2*h)
        k4_d = mult(f(y_vec_double[0] + 2*h, y_vec_double[1] + k3_d[1], y_vec_double[2] + k3_d[2]), 2*h)
        
        sum_k_d = summ(k1_d, summ(mult(k2_d, 2), summ(mult(k3_d, 2), k4_d)))
        y_vec_double = summ(y_vec_double, mult(sum_k_d, 1/6))
        
        runge_romberg_err = abs(y_vec[1] - y_vec_double[1]) / 15
        exact_err = abs(y_vec[1] - ideal_ans(y_vec[0]))
        
        print(f"In x = {y_vec[0]:.2f} | exact_error = {exact_err:.6f} | Runge-Romberg = {runge_romberg_err:.6f}")
    
    return ans


print("\n=== Runge-Kutta Method ===")
rk_result = runge_kutta(2, 7, 5, h)


# Метод Адамса 4-го порядка
def adams(x0, y0, z0, h):
    ans = [(x0, y0)]
    y_vec = [x0, y0, z0]
    f_values = []
    
    print("Starting values (Runge-Kutta):")
    for i in range(4):
        k1 = mult(f(y_vec[0], y_vec[1], y_vec[2]), h)
        k2 = mult(f(y_vec[0] + h/2, y_vec[1] + k1[1]/2, y_vec[2] + k1[2]/2), h)
        k3 = mult(f(y_vec[0] + h/2, y_vec[1] + k2[1]/2, y_vec[2] + k2[2]/2), h)
        k4 = mult(f(y_vec[0] + h, y_vec[1] + k3[1], y_vec[2] + k3[2]), h)

        sum_k = summ(k1, summ(mult(k2, 2), summ(mult(k3, 2), k4)))
        y_vec = summ(y_vec, mult(sum_k, 1/6))
        
        ans.append((y_vec[0], y_vec[1]))
        current_f = f(y_vec[0], y_vec[1], y_vec[2])
        f_values.append(current_f)
        
        exact_err = abs(y_vec[1] - ideal_ans(y_vec[0]))
        print(f"In x = {y_vec[0]:.2f} | exact_error = {exact_err:.6f}")
    
    y_vec_double = [x0, y0, z0]
    f_values_double = []
    
    for i in range(4):
        k1 = mult(f(y_vec_double[0], y_vec_double[1], y_vec_double[2]), 2*h)
        k2 = mult(f(y_vec_double[0] + h, y_vec_double[1] + k1[1]/2, y_vec_double[2] + k1[2]/2), 2*h)
        k3 = mult(f(y_vec_double[0] + h, y_vec_double[1] + k2[1]/2, y_vec_double[2] + k2[2]/2), 2*h)
        k4 = mult(f(y_vec_double[0] + 2*h, y_vec_double[1] + k3[1], y_vec_double[2] + k3[2]), 2*h)

        sum_k = summ(k1, summ(mult(k2, 2), summ(mult(k3, 2), k4)))
        y_vec_double = summ(y_vec_double, mult(sum_k, 1/6))
        
        current_f = f(y_vec_double[0], y_vec_double[1], y_vec_double[2])
        f_values_double.append(current_f)
    
    # y_vec в 2.4, y_vec_double в 2.8
    print("\nAdams-Bashforth steps:")
    for i in range(6):
        term1 = mult(f_values[-1], 55)
        term2 = mult(f_values[-2], -59)
        term3 = mult(f_values[-3], 37)
        term4 = mult(f_values[-4], -9)
        
        sum_terms = summ(term1, summ(term2, summ(term3, term4)))
        dy_new = mult(sum_terms, h / 24)
        y_vec = summ(y_vec, dy_new)
        
        ans.append((y_vec[0], y_vec[1]))
        f_new = f(y_vec[0], y_vec[1], y_vec[2])
        f_values.append(f_new)
        
        exact_err = abs(y_vec[1] - ideal_ans(y_vec[0]))
        
        runge_romberg_err = 0.0
        rr_str = "N/A"
        
        if i == 5:
            # Делаем шаг Адамса из 2.8 в 3.0
            term1_d = mult(f_values_double[-1], 55)
            term2_d = mult(f_values_double[-2], -59)
            term3_d = mult(f_values_double[-3], 37)
            term4_d = mult(f_values_double[-4], -9)
            
            sum_terms_d = summ(term1_d, summ(term2_d, summ(term3_d, term4_d)))
            dy_new_d = mult(sum_terms_d, 2*h / 24)
            y_vec_double = summ(y_vec_double, dy_new_d)
            
            runge_romberg_err = abs(y_vec[1] - y_vec_double[1]) / 15
            rr_str = f"{runge_romberg_err:.6f}"

        print(f"In x = {y_vec[0]:.2f} | exact_error = {exact_err:.6f} | Runge-Romberg = {rr_str}")
    
    return ans


print("\n=== Adams Method ===")
adams_result = adams(2, 7, 5, h)

x_euler = [p[0] for p in eiler_result]
y_euler = [p[1] for p in eiler_result]

x_rk = [p[0] for p in rk_result]
y_rk = [p[1] for p in rk_result]

x_adams = [p[0] for p in adams_result]
y_adams = [p[1] for p in adams_result]

x_exact = [2 + i*0.01 for i in range(101)]
y_exact = [ideal_ans(x) for x in x_exact]

plt.figure(figsize=(10, 6))
plt.plot(x_exact, y_exact, 'k', linewidth=3, label='Exact solution')
plt.plot(x_euler, y_euler, 'r', label='Euler', alpha=0.7)
plt.plot(x_rk, y_rk, 'b', label='Runge-Kutta', alpha=0.7)
plt.plot(x_adams, y_adams, 'g', label='Adams', alpha=0.7)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparison of numerical methods')
plt.legend()
plt.grid(True)
plt.show()