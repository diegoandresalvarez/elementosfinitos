# -*- coding: utf-8 -*-

# %% Funciones de forma del elemento rectangular lagrangiano de 16 nodos 

from mpl_toolkits.mplot3d import Axes3D  

import matplotlib.pyplot as plt
from matplotlib import cm
from sympy.polys.polyfuncs import interpolate

import numpy as np
import sympy as sp

X, Y = 0, 1

# %% Calculo las funciones de forma unidimensionales
xi, eta = sp.symbols('xi eta')
L4_xi = 4 * [None] # contenedor para las funciones de forma (en dir XI)
L4_xi[0] = sp.nsimplify(interpolate([(-1,1), (-1/3,0), (1/3,0), (1,0)], xi))
L4_xi[1] = sp.nsimplify(interpolate([(-1,0), (-1/3,1), (1/3,0), (1,0)], xi))
L4_xi[2] = sp.nsimplify(interpolate([(-1,0), (-1/3,0), (1/3,1), (1,0)], xi))
L4_xi[3] = sp.nsimplify(interpolate([(-1,0), (-1/3,0), (1/3,0), (1,1)], xi))

L4_eta = 4 * [None] # contenedor para las funciones de forma (en dir ETA)
L4_eta[0] = sp.nsimplify(interpolate([(-1,1), (-1/3,0), (1/3,0), (1,0)], eta))
L4_eta[1] = sp.nsimplify(interpolate([(-1,0), (-1/3,1), (1/3,0), (1,0)], eta))
L4_eta[2] = sp.nsimplify(interpolate([(-1,0), (-1/3,0), (1/3,1), (1,0)], eta))
L4_eta[3] = sp.nsimplify(interpolate([(-1,0), (-1/3,0), (1/3,0), (1,1)], eta))

# %% Coordenadas de los nodos y numeración local
#
#        ^ eta
#        |
#        |
# 10---9---8---7
#  |   |   |   |
# 11--16--15---6
#  |   |   |   |----> xi
# 12--13--14---5
#  |   |   |   |
#  1---2---3---4

#                  xi   eta       # nodo   
nod = np.array([[ -1  , -1   ],    #  1
                [ -1/3, -1   ],    #  2
                [  1/3, -1   ],    #  3
                [  1  , -1   ],    #  4
                [  1  , -1/3 ],    #  5
                [  1  ,  1/3 ],    #  6
                [  1  ,  1   ],    #  7
                [  1/3,  1   ],    #  8
                [ -1/3,  1   ],    #  9
                [ -1  ,  1   ],    # 10
                [ -1  ,  1/3 ],    # 11
                [ -1  , -1/3 ],    # 12
                [ -1/3, -1/3 ],    # 13
                [  1/3, -1/3 ],    # 14
                [  1/3,  1/3 ],    # 15
                [ -1/3,  1/3 ]])   # 16

nno = nod.shape[0]

# Equivalencia entre coordenada y polinomio
pos = np.empty_like(nod, dtype=int)
pos[nod==-1  ] = 0
pos[nod==-1/3] = 1
pos[nod== 1/3] = 2
pos[nod== 1  ] = 3

# %% Se calculan las funciones de forma bidimensionales
N = nno * [None]
for i in range(nno):
   N[i] = sp.simplify(L4_xi[pos[i,X]]*L4_eta[pos[i,Y]])

# Imprimo las funciones de forma
print(f'Funciones de forma lagrangianas del elemento rectangular de {nno} nodos:')
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

# %% grafico las funciones de forma
XXI  = np.linspace(-1, 1, 50)
EETA = np.linspace(-1, 1, 50)
XI, ETA = np.meshgrid(XXI, EETA)

# calculo las esferitas
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xsp = 0.025*np.cos(u)*np.sin(v)
ysp = 0.025*np.sin(u)*np.sin(v)
zsp = 0.025*np.cos(v)

for i in range(nno):
   fig = plt.figure()     # creo un lienzo
   ax = fig.gca(projection='3d')
   ax.set_aspect("equal")
   plt.xlabel(r'$\xi$', fontsize=16)   # titulo eje X
   plt.ylabel(r'$\eta$', fontsize=16)  # titulo eje Y
   plt.title(f'$N_{{{i+1}}}(\\xi,\\eta)$', fontsize=20)   

   # con este comando convierto la funcion de forma de tipo simbólico a
   # tipo función
   NN = sp.lambdify((xi, eta), N[i], 'numpy')
   surf = ax.plot_surface(XI, ETA, NN(XI,ETA), cmap=cm.viridis)
   
   # se grafican las esferitas en cada nodo
   for j in range(nno):
       ax.plot_surface(xsp+nod[j,X], ysp+nod[j,Y], zsp+(i==j), color='k')

   plt.show()

# %% bye, bye!