# -*- coding: utf-8 -*-
'''
Programa para encontrar las cargas nodales equivalentes asociados a una carga
trapezoidal:

   /|                                      ^  q(x)
   /|                             ^  |  |  |  |\
   /|                    ^  |  |  |  |  |  |  |\
   /|           ^  |  |  |  |  |  |  |  |  |  |\
   /|  ^  |  |  |  |  |  |  |  |  |  |  |  |  |\
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
   /|#########################################|\
   /|                                         |\

    |--------------------L--------------------|

Para una viga de Timoshenko, asumiendo E, I, G, A constantes
'''

# %%Importación de librerías
from sympy import symbols, integrate, solve, simplify, pretty


# %%Definción de variables
x, L, w1, w2, EI, GAast = symbols('x L w1 w2 EI GAast')
w, t, M, V = symbols('w, t, M, V')
C1, C2, C3, C4 = symbols('C1, C2, C3, C4')

q = (w2-w1)*x/L + w1  # Carga trapezoidal

# %%Solución de la ecuación diferencial de Timoshenko.
V = integrate(q, x) + C1  # Ecuaciónes calculadas.
M = integrate(V, x) + C2
t = (integrate(M, x) + C3)/EI
w = integrate(t - V/GAast, x) + C4

sol = solve([w.subs(x, 0),  # Condiciones de frontera.
             w.subs(x, L),
             t.subs(x, 0),
             t.subs(x, L)],
            [C1, C2, C3, C4])

constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

print(f'\nq = \n{pretty(q)}')
print(f'\nw = \n{pretty(simplify(w.subs(constantes)))}')
print(f'\nt = \n{pretty(simplify(t.subs(constantes)))}')
print(f'\nM = \n{pretty(simplify(M.subs(constantes)))}')
print(f'\nV = \n{pretty(simplify(V.subs(constantes)))}')

# %%Se evalúan las cargas nodales equivalentes
print(f'\nY1 = \n{pretty(simplify(-V.subs(constantes).subs(x, 0)))}')
print(f'\nY2 = \n{pretty(simplify(+V.subs(constantes).subs(x, L)))}')

# %%Y los momentos nodales equivalentes
print(f'\nM1 = \n{pretty(simplify(+M.subs(constantes).subs(x, 0)))}')
print(f'\nM2 = \n{pretty(simplify(-M.subs(constantes).subs(x, L)))}')

# Fin del programa
