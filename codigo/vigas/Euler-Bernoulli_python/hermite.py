# -*- coding: utf-8 -*-

# Interpolacion polinomica de Hermite

# FECHA         QUIEN  QUE 
# Ago 12, 2022  DAAM   El código es igual que el de MATLAB
#
# DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co

# %%Importación de librerías
import sympy as sp
from sympy.plotting import plot

# %%Definición de variables
xi, u1, u2, u1p, u2p = sp.symbols('xi u1 u2 u1p u2p')
xi1, xi2 = -1, +1

# %%Funciones de forma lagrangianas
L1 = sp.interpolate([(-1, 1), (1, 0)], xi)
L2 = sp.interpolate([(-1, 0), (1, 1)], xi)

# %%Se calculan los polinomios Pi y Qi
P1 = (1 - 2*(sp.diff(L1, xi)).subs(xi, xi1)*(xi - xi1))*L1**2
P2 = (1 - 2*(sp.diff(L2, xi)).subs(xi, xi2)*(xi - xi2))*L2**2
Q1 = (xi - xi1)*L1**2
Q2 = (xi - xi2)*L2**2

# %%Se calcula el polinomio interpolador
u = sp.simplify(P1*u1 + P2*u2 + Q1*u1p + Q2*u2p)

# %%Se calculan las funciones de forma
N1  = sp.expand(u.subs([(u1, 1), (u1p, 0), (u2, 0), (u2p, 0)])) # = P1
N1b = sp.expand(u.subs([(u1, 0), (u1p, 1), (u2, 0), (u2p, 0)])) # = Q1
N2  = sp.expand(u.subs([(u1, 0), (u1p, 0), (u2, 1), (u2p, 0)])) # = P2
N2b = sp.expand(u.subs([(u1, 0), (u1p, 0), (u2, 0), (u2p, 1)])) # = Q2

# %%Se muestran los resultados
print(f'N_1(xi)  = \n{sp.pretty(N1)}')
print(f'Nb_1(xi) = \n{sp.pretty(N1b)}')
print(f'N_2(xi)  = \n{sp.pretty(N2)}')
print(f'Nb_2(xi) = \n{sp.pretty(N2b)}')

# %%Se dibujan las funciones de forma Hermitianas
p = plot(N1, N1b, N2, N2b, (xi, -1, 1), 
         legend=True, 
         title='Funciones de forma de la viga de Euler-Bernoulli de dos nodos',
         xlabel=r'$\xi$',
         ylabel='',
         aspect_ratio=(1,1),
         show=False)
p[0].line_color = 'r'; p[0].label=r'$N_1(\xi)$'
p[1].line_color = 'b'; p[1].label=r'$\bar{N}_1(\xi)$'
p[2].line_color = 'g'; p[2].label=r'$N_2(\xi)$'
p[3].line_color = 'c'; p[3].label=r'$\bar{N}_2(\xi)$'
p.show()

# %% Bye, bye!!!