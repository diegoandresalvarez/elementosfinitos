# -*- coding: utf-8 -*-
# Programa para deducir la matriz de rigidez de un elemento de viga de
# Timoshenko a partir de la solucion de la ecuacion diferencial

# %%Importación de librerías
from sympy import symbols, integrate, pretty, simplify, zeros, solve, limit
from sympy import Matrix

# %%Definción de variables
x, L, V, M, t, w, EI, GAast = symbols('x L  V M t w EI GAast')
C1, C2, C3, C4 = symbols('C1, C2, C3, C4')  # Constantes de integración

q = 0

# %%Se calcula la matrix de rigidez
K_T = zeros(4)

V = integrate(q, x) + C1  # Ecuaciónes calculadas con constantes.
M = integrate(V, x) + C2
t = (integrate(M, x) + C3)/EI
w = integrate(t - V/GAast, x) + C4

for i in range(4):

    sol = solve([w.subs(x, 0) - int((i == 0)),  # Condiciones de frontera
                 t.subs(x, 0) - int((i == 1)),
                 w.subs(x, L) - int((i == 2)),
                 t.subs(x, L) - int((i == 3))],
                [C1, C2, C3, C4])

    constantes = [(C1, sol[C1]), (C2, sol[C2]), (C3, sol[C3]), (C4, sol[C4])]

    K_T[:, i] = [+ (V.subs(constantes)).subs(x, 0),   # Y1  se evaluan las
                 - (M.subs(constantes)).subs(x, 0),   # M1  reacciones verticales
                 - (V.subs(constantes)).subs(x, L),   # Y2  y los momentos en los
                 + (M.subs(constantes)).subs(x, L)]   # M2  apoyos

# %%Se imprime la solución
beta = (12 * EI)/(L**2 * GAast)
tmp = EI/((1 + beta)*L**3)
print(f'K_EB = EI/((1 + beta)*L**3)) * \n{pretty(simplify(K_T/tmp))}')
# Nota: se puede demostrar que:
# K22 = K44 = simplify((4+beta)*L**2)
# K24 = K42 = simplify((2-beta)*L**2)

# Se calcula la matriz de rigidez de Euler-Bernoulli
# Observe que cuando GAast -> Inf, la matriz de rigidez K se vuelve la 
# misma matriz de rigidez K de la teoria de Euler-Bernoulli:
#K_EB = limit(K_T, GAast, float('Inf'))
#print('K_EB = (EI/L^3) * \n{pretty(simplify(K_EB/(EI/L**3)))}')

# %%
#Se calculan las fuerzas nodales equivalentes:
w = symbols('w')
# q = w
q = -w*x/L

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

f = Matrix([-(V.subs(constantes)).subs(x, 0),   # Yi  se evaluan las reacciones verticales y los
            +(M.subs(constantes)).subs(x, 0),   # Mi  momentos en los apoyos y se les multiplica
            +(V.subs(constantes)).subs(x, L),   # Yj  por -1 para estimar la fuerza nodal 
            -(M.subs(constantes)).subs(x, L)])  # Mj  equivalente


beta = symbols('beta')
f = simplify(f.subs(GAast, (12 * EI)/(L**2 * beta)))
print(f'\nf = \n{pretty(f)}')

# Fin del programa
