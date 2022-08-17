# -*- coding: utf-8 -*-
'''
Programa para encontrar las cargas nodales equivalentes asociados a una carga
trapezoidal:

   /|                                      ^  q(x)
   /|                             ^  |  |  |  |/
   /|                    ^  |  |  |  |  |  |  |/
   /|           ^  |  |  |  |  |  |  |  |  |  |/
   /|  ^  |  |  |  |  |  |  |  |  |  |  |  |  |/
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/
   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/
   /|#########################################|/
   /|                                         |/

    |--------------------L--------------------|

Para una viga de Euler-Bernoulli, asumiendo E, I, A constante
'''

# FECHA         QUIEN  QUE 
# Ago 12, 2022  DAAM   El código es igual que el de MATLAB
#
# DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co

# %%Importación de librerías
import sympy as sp

# %%Definción de variables
x, L, w1, w2, EI = sp.symbols('x L w1 w2 EI')
w, t, M, V = sp.symbols('w, t, M, V')
C1, C2, C3, C4 = sp.symbols('C1, C2, C3, C4')

q = (w2-w1)*x/L + w1  # Carga trapezoidal

# %%Solución de la ecuación diferencial de Euler-Bernoulli
V =  sp.integrate(q, x) + C1  # Ecuaciónes calculadas.
M =  sp.integrate(V, x) + C2
t = (sp.integrate(M, x) + C3)/EI
w =  sp.integrate(t, x) + C4

sol = sp.solve([w.subs(x, 0),  # Condiciones de frontera.
                w.subs(x, L),
                t.subs(x, 0),
                t.subs(x, L)],
               [C1, C2, C3, C4])

constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

# se reemplazan las constantes de integración
V = V.subs(constantes)
M = M.subs(constantes)
t = t.subs(constantes)
w = w.subs(constantes)

# se imprimen los resultados
print(f'\nq = \n{sp.pretty(q)}')
print(f'\nV = \n{sp.pretty(sp.simplify(V))}')
print(f'\nM = \n{sp.pretty(sp.simplify(M))}')
print(f'\nt = \n{sp.pretty(sp.simplify(t))}')
print(f'\nw = \n{sp.pretty(sp.simplify(w))}')

# %%Se evalúan las cargas nodales equivalentes
print(f'\nY1 = \n{sp.pretty(sp.simplify(-V.subs(x, 0)))}')
print(f'\nY2 = \n{sp.pretty(sp.simplify(+V.subs(x, L)))}')

# %%Y los momentos nodales equivalentes
print(f'\nM1 = \n{sp.pretty(sp.simplify(+M.subs(x, 0)))}')
print(f'\nM2 = \n{sp.pretty(sp.simplify(-M.subs(x, L)))}')

# Fin del programa
