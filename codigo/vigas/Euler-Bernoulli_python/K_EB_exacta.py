# -*- coding: utf-8 -*-
# Programa para deducir la matriz de rigidez de un elemento de viga de
# Euler-Bernoulli a partir de la solucion de la ecuacion diferencial

# %%Importación de librerías
import sympy as sp

# %%Definción de variables
x, L, V, M, t, w, EI = sp.symbols('x L  V M t w EI')
C1, C2, C3, C4 = sp.symbols('C1, C2, C3, C4')  # Constantes de integración

q = 0

# %%Se calcula la matrix de rigidez
K = sp.zeros(4)

V =  sp.integrate(q, x) + C1  # Ecuaciónes calculadas con constantes.
M =  sp.integrate(V, x) + C2
t = (sp.integrate(M, x) + C3)/EI
w =  sp.integrate(t, x) + C4

for i in range(4):

    sol = sp.solve([w.subs(x, 0) - int(i == 0),  # Condiciones de frontera
                    t.subs(x, 0) - int(i == 1),
                    w.subs(x, L) - int(i == 2),
                    t.subs(x, L) - int(i == 3)],
                   [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K[:, i] = [+ (V.subs(constantes)).subs(x, 0),   # Y1  se evaluan las
               - (M.subs(constantes)).subs(x, 0),   # M1  reacciones verticales
               - (V.subs(constantes)).subs(x, L),   # Y2  y los momentos en los
               + (M.subs(constantes)).subs(x, L)]  # M2  apoyos

# %%Se imprime la solución
print(f'K_EB = (EI/L**3) * \n{sp.pretty(K/(EI/L**3))}')

# Fin del programa
