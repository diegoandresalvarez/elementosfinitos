# -*- coding: utf-8 -*-

# %% Deducción de las funciones de forma del elemento triangular de tres nodos

import sympy as sp
from sympy.matrices import Matrix

# %% Se definen las variables simbólicas
a1, a2, a3, x1, x2, x3, y1, y2, y3, u1, u2, u3, x, y, Area, =\
                      sp.symbols('a1 a2 a3 x1 x2 x3 y1 y2 y3 u1 u2 u3 x y Area')

# %% Resuelve el sistema de ecuaciones por a1, a2 y a3
r = sp.solve((sp.Eq(u1, a1 + a2*x1 + a3*y1),
              sp.Eq(u2, a1 + a2*x2 + a3*y2),
              sp.Eq(u3, a1 + a2*x3 + a3*y3)),  (a1, a2, a3))

# %% Muestra a1, a2 y a3
print('\na1 = '); sp.pprint(r[a1])
print('\na2 = '); sp.pprint(r[a2])
print('\na3 = '); sp.pprint(r[a3])

# %% Establece de nuevo u y se muestra
u = sp.simplify(r[a1] + r[a2]*x + r[a3]*y)
print('\nu = '); sp.pprint(u) 

# %% Se separa u en su numerador y su denominador
[num,den] = sp.fraction(u)

# %% Nosotros sabemos que el denominador es dos veces el area del triangulo.
# A continuacion se verificara:
A = sp.det(Matrix([[ 1, x1, y1 ],     # Area del triangulo con vertices
                   [ 1, x2, y2 ],     # (x1,y1), (x2,y2) y (x3,y3) numerados en el
                   [ 1, x3, y3 ]]))/2 # sentido horario de las manecillas del reloj

print(sp.simplify(sp.Eq(den, 2*A)))   # el TRUE en la respuesta confirma que el determinante es 2*Area

# %% Se vuelve a escribir u, pero con el denominador expresado como 2*Area
# ya que en el paso anterior establecimos que den es igual a 2*A
u = sp.expand(num)/(2*Area)

# %% Se factoriza u1, u2 y u3
u = sp.collect(u, (u1,u2,u3))     # Se factoriza u1, u2 y u3

# %% Se muestra # u(x,y) = N1(x,y) u1 + N2(x,y) u2 + N3(x,y) u3
print('u = '); sp.pprint(u)    # Observe aquí las funciones de forma
