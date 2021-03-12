# -*- coding: utf-8 -*-

# Programa para ilustrar que en los puntos de una cuadratura de
# Gauss-Legendre de orden n, un polinomio de grado n y otro de grado n-1
# obtenido como ajuste por mínimos cuadrados del anterior, toman el mismo
# valor.

# %%Importación de librerías

import sympy as sp
import numpy as np
from sympy.plotting import plot


# %%Interseccion en las raices del polinomio de Legendre de grado 2.

# Se definen las variables simbólicas
x, a, b = sp.symbols('x a b')

# El polinomio f y su aproximacion g
f = 1 + x + x**2
g = a + b*x

E = sp.integrate((f-g)**2, (x, -1, 1))  # Ajuste por mínimos cuadrados
r = sp.solve([sp.diff(E, a),            # diff(E,a) = 0
              sp.diff(E, b)],           # diff(E,b) = 0
             [a, b])

g = r[a] + r[b]*x         # Defino de nuevo la funcion g
inter = sp.solve(f-g, x)  # Calculo la interseccion de los polinomios

# Gráfica
p = plot(f, g, (x, -1, 1), 
         legend=True, 
         title='Intersección en las raices del polinomio de Legendre de grado 2',
         xlabel=r'$\xi$',
         ylabel='',
         aspect_ratio = (1,1),
         show=False)
p[0].line_color = 'r'; p[0].label=r'$f(\xi)$'
p[1].line_color = 'b'; p[1].label=r'$g(\xi)$'

x1, f1 = float(inter[0]), float(f.subs(x, inter[0]))
x2, f2 = float(inter[1]), float(f.subs(x, inter[1]))
p.annotations = [
{
  's'    : f'({x1:.4f}, {f1:.4f})',  # Nota en las versiones nuevas es 'text' en vez de 's'
  'xy'   : (x1 - 0.41, f1 + 0.21),
},
{
  's'    : f'({x2:.4f}, {f2:.4f})',
  'xy'   : (x2 - 0.41, f2 + 0.21),
}
]
p.show()


# %%Interseccion en las raices del polinomio de Legendre de grado 3.
import matplotlib.pyplot as plt


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
f = sp.lambdify([x], f)        # Expresión en forma de función
g = sp.lambdify([x], g)
x = np.arange(-1, 1, 0.05)  # Puntos de evaluación de las funciones.

fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(x, f(x), color='b', label=r'$f(\xi)$')
ax.plot(x, g(x), color='r', label=r'$g(\xi)$')
ax.plot(inter, f(inter), 'o', color='r', ms=10)
for i in range(3):
    xx, fx = float(inter[i]), float(f(inter[i]))
    ax.annotate(f'({xx:.4f}, {fx:.4f})', xy=(xx - 0.3, fx + 0.3))

ax.set_title('Intersección en las raices del polinomio de Legendre de grado 3')
ax.set_xlabel(r'$\xi$')
ax.legend()
ax.grid()
plt.show()

#%% Fin del programa.
