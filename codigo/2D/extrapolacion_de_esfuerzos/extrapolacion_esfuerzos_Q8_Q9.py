# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp

# %% Se define el elemento finito
EF = 9 # = {8, 9}

# %%
# Numeración local del EF serendípito rectangular de 8 nodos:
#      ^ eta
#      |
#      |
#  7---6---5
#  |   |   |
#  8---+---4----> xi
#  |   |   |
#  1---2---3

# Coordenadas de los nodos
#                   xi  eta      # nodo   
xnod8 = np.array([[ -1,  -1 ],   #  1
                  [  0,  -1 ],   #  2
                  [  1,  -1 ],   #  3
                  [  1,   0 ],   #  4
                  [  1,   1 ],   #  5
                  [  0,   1 ],   #  6
                  [ -1,   1 ],   #  7
                  [ -1,   0 ]])  #  8

# %% Numeración local del EF lagrangiano rectangular de 9 nodos:
#      ^ eta
#      |
#      |
#  7---6---5
#  |   |   |
#  8---9---4----> xi
#  |   |   |
#  1---2---3

# Coordenadas de los nodos
#                   xi  eta      # nodo   
xnod9 = np.array([[ -1,  -1 ],   #  1
                  [  0,  -1 ],   #  2
                  [  1,  -1 ],   #  3
                  [  1,   0 ],   #  4
                  [  1,   1 ],   #  5
                  [  0,   1 ],   #  6
                  [ -1,   1 ],   #  7
                  [ -1,   0 ],   #  8
                  [  0,   0 ]])  #  9

# %%
if EF == 8:
    xnod = xnod8
    mensaje = 'EF serendipito rectangular de 8 nodos'
elif EF == 9:
    xnod = xnod9
    mensaje = 'EF lagrangiano rectangular de 9 nodos'
else:
    raise Exception('Solo se permiten los EFS rectangulares de 8 y 9 nodos')

nno = xnod.shape[0]

# %% Se define la cuadratura de Gauss Legendre a utilizar
n_gl = 2
x_gl, w_gl = np.polynomial.legendre.leggauss(n_gl)

# %% número de términos del polinomio interpolador
nterm = 4 # 1   xi_gl   eta_gl   xi_gl*eta_gl

# %%  Se define la matriz A1
A1 = np.empty((nterm, nterm))

i = 0
for p in range(n_gl):
    for q in range(n_gl):
        xi_gl, eta_gl = x_gl[p], x_gl[q]
        A1[i,:] = np.array([ 1, xi_gl, eta_gl, xi_gl*eta_gl ])
        i += 1

# %% Se define la matriz A2
A2 = np.empty((nno, nterm))
for i in range(nno):
    xi_gl, eta_gl = xnod[i]
    A2[i,:] = np.array([ 1, xi_gl, eta_gl, xi_gl*eta_gl ])

# %% Se reporta la matriz
print(f'La matriz de interpolación del {mensaje} es:')
A = A2 @ np.linalg.inv(A1)    # A2/A1 de MATLAB
print(A)
