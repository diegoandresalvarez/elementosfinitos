# -*- coding: utf-8 -*-

# %% Funciones de forma del elemento rectangular serendípito de 4 y 8 nodos

import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
import sympy as sp

X, Y = 0, 1

nno = 8 #  escoja entre {4, 8}.

# %% coordenadas de los nodos y numeración local
if nno == 4:
   #      ^ eta
   #      |
   #      |
   #  4-------3
   #  |   |   |
   #  ----+---|----> xi
   #  |   |   |
   #  1-------2
   #                 xi  eta     # nodo
   nod = np.array([[ -1, -1 ],   #  1
                   [  1, -1 ],   #  2
                   [  1,  1 ],   #  3
                   [ -1,  1 ]])  #  4
elif nno == 8:
   #      ^ eta
   #      |
   #      |
   #  7---6---5
   #  |   |   |
   #  8---+---4----> xi
   #  |   |   |
   #  1---2---3
   #                 xi  eta      # nodo
   nod = np.array([[ -1, -1 ],    #  1
                   [  0, -1 ],    #  2
                   [  1, -1 ],    #  3
                   [  1,  0 ],    #  4
                   [  1,  1 ],    #  5
                   [  0,  1 ],    #  6
                   [ -1,  1 ],    #  7
                   [ -1,  0 ]])   #  8

# %% se calculan las funciones de forma bidimensionales
xxi  = nod[:, X]
eeta = nod[:, Y]
if nno == 4:
   A = np.c_[np.ones(4, dtype=int), xxi, eeta, xxi*eeta]
elif nno == 8:
   A = np.c_[np.ones(8, dtype=int), xxi, eeta,  xxi**2, xxi*eeta, eeta**2, xxi**2*eeta, xxi*eeta**2]

N = nno * [None]
xi, eta = sp.symbols('xi eta')
for i in range(nno):
   # se arma el sistema de ecuaciones
   b = np.zeros(nno);   b[i] = 1
   coef_alpha = np.linalg.solve(A, b)
   if nno == 4:
      tmp = sp.Matrix([ 1, xi, eta, xi*eta ])
   elif nno == 8:
      tmp = sp.Matrix([ 1, xi, eta, xi**2, xi*eta, eta**2, xi**2*eta, xi*eta**2 ])
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
XI, ETA = np.ogrid[-1:1:30j, -1:1:30j]

# calculo las esferitas
u, v = np.ogrid[0:2*np.pi:20j, 0:np.pi:10j]
xsp = 0.025*np.cos(u)*np.sin(v)
ysp = 0.025*np.sin(u)*np.sin(v)
zsp = 0.025*np.cos(v)

for i in range(nno):
   fig = plt.figure()     # creo un lienzo
   ax = fig.add_subplot(projection='3d')
   ax.set_xlabel(r'$\xi$', fontsize=16)   # titulo eje X
   ax.set_ylabel(r'$\eta$', fontsize=16)  # titulo eje Y
   ax.set_title(f'$N_{i+1}(\\xi,\\eta)$', fontsize=20)

   # con este comando convierto la funcion de forma de tipo simbólico a
   # tipo función
   NN = sp.lambdify((xi, eta), N[i], 'numpy')
   surf = ax.plot_surface(XI, ETA, NN(XI,ETA), cmap=cm.viridis)

   # se grafican las esferitas en cada nodo
   for j in range(nno):
       ax.plot_surface(xsp+nod[j,X], ysp+nod[j,Y], zsp+(i==j), color='k')

   plt.show()

# %% bye, bye!