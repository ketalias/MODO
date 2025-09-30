import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x, l1 = sp.symbols('x l1')


# Вирази для y та z через x, λ1
y_expr = (x + 128*l1) / (1 - 64*l1)
z_expr = -4*x / l1

# Цільова функція
F = 2*x*y_expr + 8*x*z_expr - y_expr*z_expr

a = 4096*l1**4 - 128*l1**3 + 65601*l1**2 - 2048*l1 + 16
b = -65536/5*l1**4 + 2048/5*l1**3 + 1264/5*l1**2
c = -14483456/25*l1**4 + 452608/25*l1**3 + 2864/25*l1**2

eqA = (-320*l1**3 + 5*l1**2 + 5125*l1 - 80)*x + (512*l1**3 + 632*l1**2)

eqB = a*x**2 + b*x + c


solutions = []
for guess_x in [-20, -5, -1, 0.5, 1, 5, 10, 20]:
    for guess_l in [-5, -2, -1, -0.5, 0.1, 0.5, 1, 2, 5]:
        try:
            sol = sp.nsolve([eqA, eqB], [x, l1], [guess_x, guess_l], tol=1e-14, maxsteps=100)
            sol = [float(sol[0]), float(sol[1])]
            if abs(sol[1]) > 1e-8 and abs(1 - 64*sol[1]) > 1e-8:
                if not any(abs(sol[0]-s[0])<1e-7 and abs(sol[1]-s[1])<1e-7 for s in solutions):
                    solutions.append(sol)
        except Exception:
            continue

X_vals, Y_vals, Z_vals, F_vals = [], [], [], []

print("Розв'язки з оцінкою F:")
for sol in solutions:
    xi, li = sol
    yi = float(y_expr.subs({x: xi, l1: li}))
    zi = float(z_expr.subs({x: xi, l1: li}))
    Fi = float(F.subs({x: xi, l1: li}))
    print(f"x={xi:.6f}, λ1={li:.6f}, y={yi:.6f}, z={zi:.6f}, F={Fi:.6f}")
    X_vals.append(xi)
    Y_vals.append(yi)
    Z_vals.append(zi)
    F_vals.append(Fi)

F_min, F_max = min(F_vals), max(F_vals)

# 4. Графік: точки у 3D
fig1 = plt.figure(figsize=(8,6))
ax1 = fig1.add_subplot(111, projection='3d')

for i in range(len(X_vals)):
    color = 'red'
    if abs(F_vals[i]-F_min) < 1e-6:
        color = 'green'  # мінімум
    elif abs(F_vals[i]-F_max) < 1e-6:
        color = 'blue'   # максимум
    ax1.scatter(X_vals[i], Y_vals[i], Z_vals[i], c=color, s=80)
    ax1.text(X_vals[i], Y_vals[i], Z_vals[i], f"F={F_vals[i]:.2f}", fontsize=8)

ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_title("Стаціонарні точки (зел. мін, син. макс)")

# 5. Графік: поверхня обмеження (еліпсоїд)
center_x = 8/5     
center_y = -2       
radius = 12    

u = np.linspace(0, 2*np.pi, 60)
phi = np.linspace(0, np.pi, 60)
u, phi = np.meshgrid(u, phi)

X = center_x + radius * np.cos(u) * np.sin(phi)
Y = center_y + (radius/8) * np.sin(u) * np.sin(phi)  
Z = radius * np.cos(phi)

fig2 = plt.figure(figsize=(8,6))
ax2 = fig2.add_subplot(111, projection='3d')

ax2.plot_surface(X, Y, Z, color='lightblue', alpha=0.3, edgecolor='none')

for i in range(len(X_vals)):
    color = 'red'
    if abs(F_vals[i]-F_min) < 1e-6:
        color = 'green'
    elif abs(F_vals[i]-F_max) < 1e-6:
        color = 'blue'
    ax2.scatter(X_vals[i], Y_vals[i], Z_vals[i], c=color, s=80)
    ax2.text(X_vals[i], Y_vals[i], Z_vals[i], f"P{i+1}", fontsize=8)

ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.set_title("Стаціонарні точки на поверхні обмеження")

plt.show()
