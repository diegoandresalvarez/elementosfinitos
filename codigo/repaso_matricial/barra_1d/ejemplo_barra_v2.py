# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp

# se definen algunas constantes que hacen el código más legible
NL1, NL2 = 0, 1


# %% funciones para el ensamblaje
def ensamblar_mat(K, idx, Ke):
    ''' Ensambla la matriz de rigidez local Ke en la matriz de rigidez global K

    Uso:
    ensamblar_mat(K, idx, Ke)

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
            
def ensamblar_vec(v, idx, ve):
    ''' Ensambla el vector local ve en el vector global v

    Uso:
    ensamblar_vec(v, idx, ve)

    Parametros de entrada:
    K   -> vector global
    idx -> grados de libertad donde debe proceder el ensamblaje
    Ke  -> vector local
    '''

    # se verifican las dimensiones de las matrices
    nfil, ncol = v.shape
    assert ncol == 1, "v debe ser un vector columna"

    nfil, ncol = ve.shape
    assert nfil == len(idx),\
            "ve debe ser vector columna y debe tener el mismo número de elementos que idx"

    # se procede con el ensamblaje
    for i in range(nfil):
        v[idx[i]] += ve[i]
            

# %% defino las variables
b, E, A, L, P = sp.symbols('b E A L P')         # define las variables simbólicas

long = [L,  L,  L/2]                       # longitud de la barra

# LaG: local a global: matriz que relaciona nodos locales y globales
# (se lee la barra x va del nodo i al nodo j)
LaG = np.array([[1, 3],         # fila = barra
                [2, 3],         # col1 = nodo global asociado a nodo local 1
                [3, 4]]) - 1    # col2 = nodo global asociado a nodo local 2

k = [E*A/longitud for longitud in long]    # (k minúscula) rigidez de cada barra

# %% Fuerzas nodales equivalentes para cada barra
# fe = sp.zeros((2,3));
fe = sp.Matrix([[ b*L/2, b*L/2, 0 ],
                [ b*L/2, b*L/2, 0 ]])

# %% Se prepara el vector de fuerzas nodales equivalentes global
f = sp.zeros(4,1)
f[3-1] = P/2          # se agregan las cargas puntuales
f[4-1] = P

# %% ensamblo la matriz de rigidez global (K mayúscula)
K = sp.zeros(4)       # separa memoria para matriz de rigidez global K
for e in range(3):    # para cada una de las barras e = 1, 2 y 3
   idx = LaG[e,:]     # extrae indices de los nodos globales de la barra e
   Ke  = k[e]*np.array([[1, -1],[-1, 1]])  # matriz de rigidez local
   ensamblar_mat(K, idx, Ke)               # ensambla matriz de rigidez
   ensamblar_vec(f, idx, fe[:,e])          # ensambla vector de fuerzas nodales equivalentes

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
fd  = f.extract(c,[0]);    fc  = f.extract(d,[0])
ac  = sp.Matrix([0, 0])   # SymPy los toma como vectores columna

# %% resuelvo el sistema de ecuaciones
ad = Kdd.solve(fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad - fd

# %% formo los vectores de desplazamientos (a) y fuerzas (q)
a = sp.zeros(4,1); q = sp.zeros(4,1)     # separo la memoria
for i in range(len(c)):
    a[c[i]] = ac[i]
    q[c[i]] = qd[i]

for i in range(len(d)):
    a[d[i]] = ad[i]
    # q[d[i]] = qc[i] = 0

# %% calculo las cargas axiales en cada barra
fax = sp.zeros(3,1)
for e in range(3):
   fax[e] = k[e]*(a[LaG[e,NL2]] - a[LaG[e,NL1]])

# %% se simplifican las expresiones
a   = sp.simplify(a)
q   = sp.simplify(q)
fax = sp.simplify(fax)

# %% imprimo los resultados
print('\n\na = \n');   sp.pprint(a)
print('\n\nq = \n');   sp.pprint(q)
print('\n\nfax = \n'); sp.pprint(fax)
