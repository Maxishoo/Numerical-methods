import matplotlib.pyplot as plt

def f(x):
    return x**2 / (x**4 + 256)

X0 = 0 
Xk = 2
h1 = 0.5
h2 = 0.25

# Вычисление для h1 = 0.5
n1 = int((Xk - X0) / h1)
X1 = [X0 + i*h1 for i in range(n1+1)]

# Метод прямоугольников (средних) для h1
integral_rec_h1 = sum([h1 * f((X1[i-1] + X1[i])/2) for i in range(1, n1+1)])

# Метод трапеций для h1
integral_trap_h1 = sum([h1 * (f(X1[i-1]) + f(X1[i])) / 2 for i in range(1, n1+1)])

# Метод Симпсона для h1 (требует чётного числа интервалов)
integral_sim_h1 = 0
for i in range(1, n1, 2):
    integral_sim_h1 += (h1/3) * (f(X1[i-1]) + 4*f(X1[i]) + f(X1[i+1]))

# Вычисление для h2 = 0.25
n2 = int((Xk - X0) / h2)
X2 = [X0 + i*h2 for i in range(n2+1)]

# Метод прямоугольников (средних) для h2
integral_rec_h2 = sum([h2 * f((X2[i-1] + X2[i])/2) for i in range(1, n2+1)])

# Метод трапеций для h2
integral_trap_h2 = sum([h2 * (f(X2[i-1]) + f(X2[i])) / 2 for i in range(1, n2+1)])

# Метод Симпсона для h2
integral_sim_h2 = 0
for i in range(1, n2, 2):
    integral_sim_h2 += (h2/3) * (f(X2[i-1]) + 4*f(X2[i]) + f(X2[i+1]))

# Точное значение (вычисленное аналитически или очень маленьким шагом)
# Для проверки используем метод Симпсона с очень малым шагом
h_exact = 0.0001
n_exact = int((Xk - X0) / h_exact)
X_exact = [X0 + i*h_exact for i in range(n_exact+1)]
exact_value = 0
for i in range(1, n_exact, 2):
    exact_value += (h_exact/3) * (f(X_exact[i-1]) + 4*f(X_exact[i]) + f(X_exact[i+1]))

# Метод Рунге-Ромберга-Ричардсона
# k = h1/h2 = 0.5/0.25 = 2
k = 2

# Прямоугольники: p = 2 (метод средних прямоугольников)
F_rec = integral_rec_h2 + (integral_rec_h2 - integral_rec_h1) / (k**2 - 1)

# Трапеции: p = 2
F_trap = integral_trap_h2 + (integral_trap_h2 - integral_trap_h1) / (k**2 - 1)

# Симпсон: p = 4
F_sim = integral_sim_h2 + (integral_sim_h2 - integral_sim_h1) / (k**4 - 1)

# Вычисление погрешностей
error_rec_h1 = abs(integral_rec_h1 - exact_value)
error_rec_h2 = abs(integral_rec_h2 - exact_value)
error_rec_RR = abs(F_rec - exact_value)

error_trap_h1 = abs(integral_trap_h1 - exact_value)
error_trap_h2 = abs(integral_trap_h2 - exact_value)
error_trap_RR = abs(F_trap - exact_value)

error_sim_h1 = abs(integral_sim_h1 - exact_value)
error_sim_h2 = abs(integral_sim_h2 - exact_value)
error_sim_RR = abs(F_sim - exact_value)

# Красивая таблица в консоли
print("=" * 95)
print(" " * 35 + "СРАВНЕНИЕ ЧИСЛЕННЫХ МЕТОДОВ ИНТЕГРИРОВАНИЯ")
print("=" * 95)
print(f"{'Метод':<20} {'h = 0.5':>20} {'h = 0.25':>20} {'Рунге-Ромберг':>20} {'Точное':>15}")
print("-" * 95)
print(f"{'Прямоугольники':<20} {integral_rec_h1:>20.10f} {integral_rec_h2:>20.10f} {F_rec:>20.10f} {exact_value:>15.10f}")
print(f"{'Трапеции':<20} {integral_trap_h1:>20.10f} {integral_trap_h2:>20.10f} {F_trap:>20.10f} {exact_value:>15.10f}")
print(f"{'Симпсон':<20} {integral_sim_h1:>20.10f} {integral_sim_h2:>20.10f} {F_sim:>20.10f} {exact_value:>15.10f}")
print("=" * 95)

print("\n" + "=" * 95)
print(" " * 35 + "ПОГРЕШНОСТИ (абсолютные)")
print("=" * 95)
print(f"{'Метод':<20} {'|Ошибка| при h=0.5':>25} {'|Ошибка| при h=0.25':>25} {'|Ошибка| Рунге-Ромберг':>25}")
print("-" * 95)
print(f"{'Прямоугольники':<20} {error_rec_h1:>25.10e} {error_rec_h2:>25.10e} {error_rec_RR:>25.10e}")
print(f"{'Трапеции':<20} {error_trap_h1:>25.10e} {error_trap_h2:>25.10e} {error_trap_RR:>25.10e}")
print(f"{'Симпсон':<20} {error_sim_h1:>25.10e} {error_sim_h2:>25.10e} {error_sim_RR:>25.10e}")
print("=" * 95)

# Визуализация
plt.figure(figsize=(10, 6))
x_plot = [i/1000 for i in range(2001)]
y_plot = [f(x) for x in x_plot]
plt.plot(x_plot, y_plot, 'b-', linewidth=2, label='f(x) = x²/(x⁴+256)')
plt.fill_between(x_plot, y_plot, alpha=0.3)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('График подынтегральной функции')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()