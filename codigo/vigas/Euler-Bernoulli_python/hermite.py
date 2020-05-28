# -*- coding: utf-8 -*-

# Interpolacion polinomica de Hermite.


# %%Importación de librerías

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# %%Definición de variables

xi, u1, u2, u1p, u2p = sp.symbols('xi u1 u2 u1p u2p')
x1 = -1
x2 =  1

# %%Funciones de forma lagrangianas

L1 = sp.interpolate([(-1, 1), (1, 0)], xi)
L2 = sp.interpolate([(-1, 0), (1, 1)], xi)

# %%Se calculan los polinomios Pi y

P1 = (1 - 2*(sp.diff(L1, xi)).subs(xi, x1)*(xi-x1))*L1**2
P2 = (1 - 2*(sp.diff(L2, xi)).subs(xi, x2)*(xi-x2))*L2**2
Q1 = (xi-x1)*L1**2
Q2 = (xi-x2)*L2**2

# %%Se calcula el polinomio interpolador

u = sp.simplify(P1*u1 + P2*u2 + Q1*u1p + Q2*u2p)

# %%Se calculan las funciones de forma
N1  = sp.expand(u.subs([(u1, 1), (u1p, 0), (u2, 0), (u2p, 0)]))
N1b = sp.expand(u.subs([(u1, 0), (u1p, 1), (u2, 0), (u2p, 0)]))
N2  = sp.expand(u.subs([(u1, 0), (u1p, 0), (u2, 1), (u2p, 0)]))
N2b = sp.expand(u.subs([(u1, 0), (u1p, 0), (u2, 0), (u2p, 1)]))

# %%Se muestran los resultados

print(f'N_1(xi)  = \n{sp.pretty(sp.pretty(N1))}')
print(f'Nb_1(xi) = \n{sp.pretty(N1b)}')
print(f'N_2(xi)  = \n{sp.pretty(N2)}')
print(f'Nb_2(xi) = \n{sp.pretty(N2b)}')

# %%Se dibujan las funciones de forma Hermitianas

N1p  = sp.lambdify([xi], N1)
N1bp = sp.lambdify([xi], N1b)
N2p  = sp.lambdify([xi], N2)
N2bp = sp.lambdify([xi], N2b)
xical = np.arange(-1, 1, 0.05)  # Puntos de evaluación de las funciones.

fig, ax = plt.subplots()
ax.plot(xical, N1p(xical), linewidth=2, label=r'$N_1(\xi)$')
ax.plot(xical, N1bp(xical), linewidth=2, label=r'$\bar{N}_1(\xi)$')
ax.plot(xical, N2p(xical), linewidth=2, label=r'$N_2(\xi)$')
ax.plot(xical, N2bp(xical), linewidth=2, label=r'$\bar{N}_2(\xi)$')

ax.axis('equal')
ax.set_title('Funciones de forma de la viga de Euler-Bernoulli de dos nodos')
ax.set_xlabel(r'$\xi$')

ax.legend()
ax.grid()
plt.show()

# Fin del programa.
