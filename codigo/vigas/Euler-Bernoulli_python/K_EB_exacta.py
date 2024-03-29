# -*- coding: utf-8 -*-
# Programa para deducir la matriz de rigidez de un elemento de viga de
# Euler-Bernoulli a partir de la solucion de la ecuacion diferencial

# FECHA         QUIEN  QUE 
# Ago 11, 2022  DAAM   El código es igual que el de MATLAB
#
# DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co


# %%Importación de librerías
import sympy as sp

# %% Definición de variables
x, L, V, M, t, w, EI = sp.symbols('x L  V M t w EI')
C1, C2, C3, C4 = sp.symbols('C1, C2, C3, C4')  # Constantes de integración

q = 0

# %% Se calcula la matrix de rigidez
K = sp.zeros(4)

V =  sp.integrate(q, x) + C1         # se integran las ec. diferenciales
M =  sp.integrate(V, x) + C2         # de la viga de Euler-Bernoulli
t = (sp.integrate(M, x) + C3)/EI
w =  sp.integrate(t, x) + C4
for i in range(4):
    sol = sp.solve([sp.Eq(w.subs(x, 0), int(i == 0)),  # Condiciones de frontera
                    sp.Eq(t.subs(x, 0), int(i == 1)),
                    sp.Eq(w.subs(x, L), int(i == 2)),
                    sp.Eq(t.subs(x, L), int(i == 3))],
                   [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K[:, i] = [+V.subs(constantes).subs(x, 0),   # Y1  se evaluan las
               -M.subs(constantes).subs(x, 0),   # M1  reacciones verticales
               -V.subs(constantes).subs(x, L),   # Y2  y los momentos en los
               +M.subs(constantes).subs(x, L)]   # M2  apoyos

# %% Se imprime la solución
print(f'K_EB = (EI/L**3) * \n{sp.pretty(K/(EI/L**3))}')

# %% Bye, bye!
