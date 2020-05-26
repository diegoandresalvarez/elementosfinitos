# -*- coding: utf-8 -*-
'''
Programa para deducir la matriz de rigidez de un elemento de viga de
Euler-Bernoulli a partir de la solucion de la ecuacion diferencial
que tiene una rótula a la derecha.
'''

# %%Importación de librerías
import sympy as sp

# %%Definción de variables
x, L, V, M, t, w, EI = sp.symbols('x L  V M t w EI')
C1, C2, C3, C4 = sp.symbols('C1, C2, C3, C4')  # Constantes de integración

q = 0

# %%Se calcula la matriz de rigidez para la rotula en el lado derecho
K = sp.zeros(3)

V = sp.integrate(q, x) + C1  # Ecuaciónes calculadas con constantes.
M = sp.integrate(V, x) + C2
t = (sp.integrate(M, x) + C3)/EI
w = sp.integrate(t, x) + C4

for i in range(3):

    sol = sp.solve([w.subs(x, 0) - int((i == 0)),  # Condiciones de frontera
                    t.subs(x, 0) - int((i == 1)),
                    w.subs(x, L) - int((i == 2)),
                    M.subs(x, L)],
                   [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K[:, i] = [+ (V.subs(constantes)).subs(x, 0),  # Y1 se evaluan las
               - (M.subs(constantes)).subs(x, 0),  # M1 reacciones vert y los
               - (V.subs(constantes)).subs(x, L)]  # Y2 momentos en los apoyos

# Se imprime la solución
print('La matriz de rigidez para una rótula en el lado derecho es:')
print(f'{sp.pretty(K)}')

# %%Se calcula la matriz de rigidez para la rotula en el lado izquierdo
K = sp.zeros(3)

V = sp.integrate(q, x) + C1  # Ecuaciónes calculadas con constantes.
M = sp.integrate(V, x) + C2
t = (sp.integrate(M, x) + C3)/EI
w = sp.integrate(t, x) + C4

for i in range(3):

    sol = sp.solve([w.subs(x, 0) - int((i == 0)),  # Condiciones de frontera
                    M.subs(x, 0),
                    w.subs(x, L) - int((i == 1)),
                    t.subs(x, L) - int((i == 2))],
                   [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K[:, i] = [+ (V.subs(constantes)).subs(x, 0),  # Y1 se evaluan las
               - (V.subs(constantes)).subs(x, L),  # Y2 reacciones vert y los
               + (M.subs(constantes)).subs(x, L)]  # M2 momentos en los apoyos

# Se imprime la solución
print('\nLa matriz de rigidez para una rótula en el lado izquierdo es:')
print(f'{sp.pretty(K)}')
