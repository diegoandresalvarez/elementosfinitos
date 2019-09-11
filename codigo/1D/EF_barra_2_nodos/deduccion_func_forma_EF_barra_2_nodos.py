# Deducción de las funciones de forma del elemento de barra de 2 nodos

import sympy as sp
from sympy.matrices import Matrix

# Defino las variables simbolicas
u1, u2, x, x1, x2, a1, a0 = sp.symbols('u1 u2 x x1 x2 a1 a0')

r = sp.solve((sp.Eq(u1, a1*x1 + a0),
              sp.Eq(u2, a1*x2 + a0)), (a0, a1))

print('a0 = '); sp.pprint(r[a0])
print('a1 = '); sp.pprint(r[a1])

u = sp.expand(r[a1]*x + r[a0]) # Se define ahora u(x) ya que conocemos a1 y a0
u = sp.collect(u, (u1,u2))     # Se factoriza u2
print('u = '); sp.pprint(u)    # Observe aquí las funciones de forma
