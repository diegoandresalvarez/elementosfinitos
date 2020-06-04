# -*- coding: utf-8 -*-

# %% Funciones de forma del elemento triangular de 10 nodos

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
import sympy as sp

X, Y = 1, 2
L1, L2, L3 = 0, 1, 2

nno = 10

# %% coordenadas de los nodos y numeración local
if nno == 10:
   #  ^ eta
   #  |
   #  |
   #  3
   #  | \
   #  8---7
   #  |   | \
   #  9--10---6
   #  |   |   | \  
   #  1---4---5---2----> xi
   #
   #                 xi  eta      do      nodo
   nod = np.array([[1  ,  0  ,  0   ],  #  1
                   [0  ,  1  ,  0   ],  #  2
                   [0  ,  0  ,  1   ],  #  3
                   [2/3,  1/3,  0   ],  #  4
                   [1/3,  2/3,  0   ],  #  5
                   [0  ,  2/3,  1/3 ],  #  6
                   [0  ,  1/3,  2/3 ],  #  7
                   [1/3,  0  ,  2/3 ],  #  8
                   [2/3,  0  ,  1/3 ],  #  9
                   [1/3,  1/3,  1/3 ]]) # 10

# %% se calculan las funciones de forma bidimensionales
xxi  = nod[:, L2]
eeta = nod[:, L3]
if nno == 10:
   A = np.c_[np.ones(nno, dtype=int), xxi, eeta,  xxi**2, xxi*eeta, eeta**2, xxi**3, xxi**2*eeta, xxi*eeta**2, eeta**3]

N = nno * [None]
xi, eta = sp.symbols('xi eta')
for i in range(nno):
   # se arma el sistema de ecuaciones
   b = np.zeros(nno);   b[i] = 1
   coef_alpha = np.linalg.solve(A, b)
   coef_alpha[np.abs(coef_alpha) < 1e-10] = 0
   if nno == 10:
      tmp = sp.Matrix([ 1, xi, eta, xi**2, xi*eta, eta**2, xi**3, xi**2*eta, xi*eta**2, eta**3 ])
   N[i] = sp.nsimplify(tmp.dot(coef_alpha))

# %% imprimo las funciones de forma
print(f'Funciones de forma serendípitas del elemento rectangular de {nno} nodos:')
for i in range(nno):
   print(f'N{i+1}(xi,eta) = {N[i]}')

# se calculan las derivadas de las funciones de forma con respecto a xi y
# con respecto a eta y se imprimen (para referencias posteriores):
print('\nDerivadas con respecto a xi:')
for i in range(nno):
   print(f'dN{i+1}(xi,eta)_dxi = {sp.diff(N[i], xi)}')

print('\nDerivadas con respecto a eta:')
for i in range(nno):
   print(f'dN{i+1}(xi,eta)_deta = {sp.diff(N[i], eta)}')

print()

# %% Se verifica la condición de cuerpo rígido:
print('\nSe verifica la condición de cuerpo rígido: sum(N) == 1:')
print(sp.simplify(sum(N)) == 1)

# %% grafico las funciones de forma
XI, ETA = np.mgrid[0:1:30j, 0:1:30j]
L1 = 1 - XI - ETA

# calculo las esferitas
u, v = np.ogrid[0:2*np.pi:20j, 0:np.pi:10j]
xsp = 0.025*np.cos(u)*np.sin(v)
ysp = 0.025*np.sin(u)*np.sin(v)
zsp = 0.025*np.cos(v)

for i in range(nno):
   fig = plt.figure()     # creo un lienzo
   ax = fig.gca(projection='3d')
   #ax.set_aspect("equal")
   plt.xlabel(r'$\xi$', fontsize=16)   # titulo eje X
   plt.ylabel(r'$\eta$', fontsize=16)  # titulo eje Y
   plt.title(f'$N_{{{i+1}}}(\\xi,\\eta)$', fontsize=20)

   # con este comando convierto la funcion de forma de tipo simbólico a
   # tipo función
   NN = sp.lambdify((xi, eta), N[i], 'numpy')
   FF = NN(XI,ETA)
   FF[L1<0] = np.NaN
   surf = ax.plot_surface(XI, ETA, FF, cmap=cm.viridis)

   # se grafican las esferitas en cada nodo
   for j in range(nno):
       ax.plot_surface(xsp+nod[j,X], ysp+nod[j,Y], zsp+(i==j), color='k')

   plt.show()

# %% bye, bye!
