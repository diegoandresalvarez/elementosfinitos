# -*- coding: utf-8 -*-
# Imposición de las deformaciones angulares gxz
# Se realiza para el elemento finito de Timoshenko de dos nodos

# %%Importación  de librerías
from sympy import Matrix, symbols, integrate, pretty, nsimplify

# %%Definción de variables
G, Aast, L, xi, w1, w2, t1, t2 = symbols('G Aast L xi w1 w2 t1 t2')

# La matriz Bs "normal"
Bs = lambda xi: Matrix([-1/L, xi/2 - 1/2, 1/L, - xi/2 - 1/2]).T

# Las funciones de interpolacion de las deformaciones gxz
Ng = 1  # se evaluara en el centro

# %%Cálculos
# Se calcula la matriz Bs_sustitutiva
Bs_sustitutiva = Ng*Bs(0)

# y la correspondiente matriz de rigidez
Ks = integrate(Bs_sustitutiva.T*G*Aast*Bs_sustitutiva*L/2, (xi, -1, 1))

print(f'Ks = (G*Aast/L) * \n{pretty(nsimplify(Ks/(G*Aast/L)))}')

# Este Ks coincide con aquel integrado con un punto de GL.

# Fin del programa
