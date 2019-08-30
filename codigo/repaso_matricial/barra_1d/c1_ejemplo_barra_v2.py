# -*- coding: utf-8 -*-

# Ejemplo 1.1. Oñate

import numpy as np
import sympy as sp

def ensamblar(K, idx, Ke):
    ''' Ensambla la matriz de rigidez local Ke en la matriz de rigidez global K

    Uso:
    ensamblar(K, idx, Ke)

    Parametros de entrada:
    K   -> matriz de rigidez global
    idx -> grados de libertad donde debe proceder el ensamblaje
    Ke  -> matriz de rigidez local
    '''

    # se verifican las dimensiones de las matrices
    nfil, ncol = K.shape
    assert nfil == ncol, "La matriz de rigidez K debe ser cuadrada"

    nfil, ncol = Ke.shape
    assert nfil == ncol == len(idx),\
            "Ke debe ser cuadrada y debe tener el mismo número de filas que idx"

    # se procede con el ensamblaje
    for i in range(nfil):
        for j in range(ncol):
            K[idx[i], idx[j]] += Ke[i,j]

# %% defino las variables
E, A, L, P = sp.symbols('E A L P')         # define las variables simbólicas

long = [L,  L,  L/2]                       # longitud de la barra

# LaG: local a global: matriz que relaciona nodos locales y globales
# (se lee la barra x va del nodo i al nodo j)
LaG = np.array([[1, 3],         # fila = barra
                [2, 3],         # col1 = nodo global asociado a nodo local 1
                [3, 4]]) - 1    # col2 = nodo global asociado a nodo local 2

k = [E*A/longitud for longitud in long]    # (k minúscula) rigidez de cada barra

# %% ensamblo la matriz de rigidez global (K mayúscula)
K = sp.zeros(4)       # separa memoria para matriz de rigidez global K
for e in range(3):    # para cada una de las barras e = 1, 2 y 3
   idx = LaG[e,:]     # extrae indices de los nodos globales de la barra e
   Ke  = k[e]*np.array([[1, -1],[-1, 1]])  # matriz de rigidez local
   ensamblar(K, idx, Ke)                   # ensamblaje matricial

print("Imprimamos la matriz de rigidez =")   
sp.pprint(K)

# %% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = np.array([1, 2]) - 1;    d = np.setdiff1d(np.arange(4), c)
c = sp.Matrix(c);            d = sp.Matrix(d)   # se convierten de Numpy a SymPy

# %% extraigo las submatrices y especifico las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |     recuerde siempre que qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |     en este caso en particular fd=0
Kcc = K.extract(c,c);      Kcd = K.extract(c,d)
Kdc = K.extract(d,c);      Kdd = K.extract(d,d)

ac = sp.Matrix([0, 0])   # SymPy los toma como vectores columna
fc = sp.Matrix([0, P])

# %% resuelvo el sistema de ecuaciones
ad = Kdd.solve(fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad

# %% formo los vectores de desplazamientos (a) y fuerzas (q)
a = sp.zeros(4,1); q = sp.zeros(4,1)     # separo la memoria
for i in range(len(c)):
    a[c[i]] = ac[i]
    q[c[i]] = qd[i]

for i in range(len(d)):
    a[d[i]] = ad[i]
    # q[d[i]] = qc[i] = 0

# %% calculo las cargas axiales en cada barra
N = sp.zeros(3,1)
for e in range(3):
   N[e] = k[e]*(a[LaG[e,2-1]] - a[LaG[e,1-1]])

# %% imprimo los resultados
print('\n\na = \n'); sp.pprint(a)
print('\n\nq = \n'); sp.pprint(q)
print('\n\nN = \n'); sp.pprint(N)
