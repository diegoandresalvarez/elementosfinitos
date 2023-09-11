# -*- coding: utf-8 -*-

# %% Programa para obtener K y f correspondientes al elemento isoparamétrico 
# lagrangiano cuadrático de tres nodos unidimensional 

import numpy as np
import sympy as sp
from sympy.polys.polyfuncs import interpolate
from sympy.matrices import Matrix
sp.init_printing(pretty_print=True)

# %% Definición de variables
xi, x1, x2, x3, L, E, A, b = sp.symbols('xi x1 x2 x3 L E A b')

# %% Defino las posiciones de los nodos
x3 = x1+L
x2 = (x1+x3)/2

# %% Funciones de forma lagrangianas
N1 = interpolate([(-1,1), (0,0), (1,0)], xi)  # = xi*(xi-1)/2
N2 = interpolate([(-1,0), (0,1), (1,0)], xi)  # = (1+xi)*(1-xi)
N3 = interpolate([(-1,0), (0,0), (1,1)], xi)  # = xi*(xi+1)/2

# %% Interpolación de la geometría y sus derivadas
x      = sp.simplify(N1*x1 + N2*x2 + N3*x3)     
dx_dxi = sp.diff(x, xi)

# NOTA: el siguiente comando solo se puede realizar si x(xi) es una función
# continua e inyectiva. Este es nuestro caso:
# https://en.wikipedia.org/wiki/Inverse_functions_and_differentiation
dxi_dx = 1/dx_dxi
# recuerde que se debe garantizar que dx_dxi>0 y dxi_dx>0

# %% Definición de la matriz de forma N y matriz de deformación del elemento B
N = Matrix([[N1, N2, N3]]) 
# B = sp.simplify([sp.diff(N1,xi) sp.diff(N2,xi) sp.diff(N3,xi)])
B = sp.diff(N,xi)*dxi_dx

# %% "matriz constitutiva"
D = E*A

# %% Calculo la matriz de rigidez del elemento
fc = (A*E)/(3*L) # factor común
K = sp.simplify(sp.integrate(B.T*D*B*dx_dxi, (xi, -1, 1)))
K = sp.MatMul(fc, K/fc, evaluate = False)
print('K = '); sp.pprint(K, num_columns=150); print()

# %% Calculo la matriz de fuerzas nodales equivalentes del elemento
fc = (b*L)/6 # factor común
f = sp.simplify(sp.integrate(N.T*b*dx_dxi, (xi, -1, 1)))
f = sp.MatMul(fc, f/fc, evaluate = False)
print('f = '); sp.pprint(f)
