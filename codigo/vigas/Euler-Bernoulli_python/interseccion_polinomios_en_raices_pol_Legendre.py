# -*- coding: utf-8 -*-

# Programa para ilustrar que en los puntos de una cuadratura de
# Gauss-Legendre de orden n, un polinomio de grado n y otro de grado n-1
# obtenido como ajuste por minimos cuadrados del anterior, toman el mismo
# valor.

# %%Importación de librerías

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# %%Interseccion en las raices del polinomio de Legendre de grado 2.

# Se definen las variables simbólicas
x, a, b = sp.symbols('x a b')

# El polinomio f y su aproximacion g
f = 1 + x + x**2
g = a + b*x

E = sp.integrate((f-g)**2, (x, -1, 1))  # El ajuste por minimos cuadrados
r = sp.solve([sp.diff(E, a),            # diff(E,a) = 0
              sp.diff(E, b)],           # diff(E,b) = 0
             [a, b])

g = r[a] + r[b]*x         # Defino de nuevo la funcion g
inter = np.array(sp.solve(f-g, x))  # Cálculo la interseccion de los polinomios

# Gráfica
fcal = sp.lambdify([x], f)        # Expresión en forma de función
gcal = sp.lambdify([x], g)
xcal = np.arange(-1, 1, 0.05)  # Puntos de evaluación de las funciones.

fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(xcal, fcal(xcal), color='b')
ax.plot(xcal, gcal(xcal), color='r')
ax.plot(inter, fcal(inter), 'o', color='r', ms=10)
for i in range(2):
    ax.annotate(f'({round(inter[i].evalf(),2)}, '
                + f'{round(fcal(inter[i].evalf()),2)})',
                xy=(round(inter[i].evalf(), 2) - 0.41,
                    round(fcal(inter[i].evalf()), 2) + 0.21),
                fontsize=13)

ax.set_title('Intersección en las raices del polinomio de Legendre de grado 2')
ax.grid()

# %%Interseccion en las raices del polinomio de Legendre de grado 3.

# Se definen las variables simbólicas
x, a, b, c = sp.symbols('x a b c')

# El polinomio f y su aproximacion g
f = 1 + x + x**2 + x**3
g = a + b*x + c*x**2

E = sp.integrate((f-g)**2, (x, -1, 1))  # El ajuste por minimos cuadrados
r = sp.solve([sp.diff(E, a),            # diff(E,a) = 0
              sp.diff(E, b),            # diff(E,b) = 0
              sp.diff(E, c)],           # diff(E,c) = 0
             [a, b, c])

g = r[a] + r[b]*x + r[c]*x**2  # Defino de nuevo la funcion g
inter = np.array(sp.solve(f-g, x))  # Cálculo la interseccion de los polinomios

# Gráfica
fcal = sp.lambdify([x], f)        # Expresión en forma de función
gcal = sp.lambdify([x], g)
xcal = np.arange(-1, 1, 0.05)  # Puntos de evaluación de las funciones.

fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(xcal, fcal(xcal), color='b')
ax.plot(xcal, gcal(xcal), color='r')
ax.plot(inter, fcal(inter), 'o', color='r', ms=10)
for i in range(3):
    ax.annotate(f'({round(inter[i].evalf(),2)}, '
                + f'{round(fcal(inter[i].evalf()),2)})',
                xy=(round(inter[i].evalf(), 2) - 0.3,
                    round(fcal(inter[i].evalf()), 2) + 0.3),
                fontsize=13)

ax.set_title('Intersección en las raices del polinomio de Legendre de grado 3')
ax.grid()
plt.show()

# Fin del programa.
