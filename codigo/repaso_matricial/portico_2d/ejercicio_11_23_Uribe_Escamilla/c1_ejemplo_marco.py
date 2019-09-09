# -*- coding: utf-8 -*-

import numpy as np

# %% Ejemplo 11.23 Uribe Escamilla

# %% Unidades en toneladas y metros

# %% define sin() y cos(), pero con un argumento en grados (como en MATLAB)
np.cosd = lambda x : np.cos(np.deg2rad(x))
np.sind = lambda x : np.sin(np.deg2rad(x))

# se definen algunas constantes que hacen el código más legible
NL1, NL2 = 0 ,1
X,   Y   = 0, 1

# %% defino las variables
Aviga = 0.30*0.35;        Acol  = 0.30*0.30        # m^2    área
Iviga = 0.30*0.35**3/12;  Icol  = 0.30*0.30**3/12  # m^4    inercia_y

# barra   1                 2                 3
A     = [ Aviga,            Acol,             Acol             ] # áreas
I     = [ Iviga,            Icol,             Icol             ] # inercias_y
long  = [ np.hypot(4,2),    5,                np.hypot(2,6)    ] # longitud barra (m)
theta = [ np.arctan2(2,4),  np.arctan2(4,3),  np.arctan2(-6,2) ] # ángulo inclinación (rad)
theta = np.degrees(theta) # rad -> grados

# LaG: local a global: matriz que relaciona nodos locales y globales
# (se lee la barra x va del nodo i al nodo j)
LaG = np.array([[1, 2],       # fila = barra
                [4, 1],       # col1 = nodo global asociado a nodo local 1
                [2, 3]]) - 1  # col2 = nodo global asociado a nodo local 2                

# gdl: grados de libertad
ngdl = 12 # numero de grados de libertad
gdl = np.array([[ 4,  5,  6],      # fila = nodo
                [ 7,  8,  9],      # col1 = gdl en dirección x
                [10, 11, 12],      # col2 = gdl en dirección y
                [ 1,  2,  3]]) - 1 # col3 = gdl en dirección angular antihoraria

E = 190*10000 # ton/m^2  modulo de elasticidad

nb = LaG.shape[0] # número de barras (número de filas de LaG)
nn = gdl.shape[0] # número de nodos  (número de filas de gdl)

# %% fuerzas nodales equivalentes para las diferentes barras
# (en este ejemplo las fuerzas nodales equivalentes estas siendo 
# especificadas con respecto al sistema de coordenadas globales)
fe = nb*[None] # = cell(nb,1) de MATLAB
#          fxi    fyi     mi      fxj    fyj     mj
#          ton    ton     ton-m   ton    ton     ton-m
fe[1-1] = [0,    -5.60,  -3.733,  0,    -5.60,   +3.733 ] # OJO con los signos
fe[2-1] = [0,     0,      0,      0,     0,       0     ] # mirar pag 613
fe[3-1] = [0,     0,      0,      0,     0,       0     ]
fe = [np.array(v) for v in fe]  # convierto todas esas listas a arrays de NumPy

# %% separo la memoria
K   = np.zeros((ngdl, ngdl))  # matriz de rigidez global
f   = np.zeros(ngdl)          # vector de fuerzas nodales equivalentes global
Ke  = nb*[None]               # matriz de rigidez local en coordenadas globales
T   = nb*[None]               # matriz de transfor(mación de coordenadas
idx = np.zeros((nb, 6), dtype=int) # almacena los 6 gdls de las barras

# %% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e in range(nb): # para cada barra
    # saco los 6 gdls de la barra e
    idx[e,:] = np.r_[gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:]]
    
    # matriz de transformaciÓn de coordenadas para la barra e
    c = np.cosd(theta[e]); s = np.sind(theta[e]);  # sin y cos de la inclinación
    T[e] = np.array([[ c,  s,  0,  0,  0,  0],
                     [-s,  c,  0,  0,  0,  0],
                     [ 0,  0,  1,  0,  0,  0],
                     [ 0,  0,  0,  c,  s,  0],
                     [ 0,  0,  0, -s,  c,  0],
                     [ 0,  0,  0,  0,  0,  1]])
          
    # matriz de rigidez local expresada en el sistema de coordenadas locales
    # para la barra e
    AE = A[e]*E;    EI = E*I[e];    L=long[e]; L2=long[e]**2; L3=long[e]**3
    Kloc = np.array(
       [[ AE/L,    0,         0,        -AE/L,    0,         0       ],
        [  0,     12*EI/L3,   6*EI/L2,   0,     -12*EI/L3,   6*EI/L2 ],
        [  0,      6*EI/L2,   4*EI/L,    0,      -6*EI/L2,   2*EI/L  ],
        [-AE/L,    0,         0,         AE/L,    0,         0       ],
        [  0,    -12*EI/L3,  -6*EI/L2,   0,      12*EI/L3,  -6*EI/L2 ],
        [  0,      6*EI/L2,   2*EI/L,    0,      -6*EI/L2,   4*EI/L  ]])
 
    # matriz de rigidez local en coordenadas globales
    Ke[e] = T[e].T@Kloc@T[e]         
 
    K[np.ix_(idx[e,:],idx[e,:])] += Ke[e] # ensamblo Ke en K global
    f[idx[e,:]]                  += fe[e] # ensamblo fe en f global

# %% localizo la carga puntual de 1.5 ton en el gdl 4
f[4-1] += 1.5

# %% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = np.array([1, 2, 3, 10, 11, 12]) - 1;
d = np.setdiff1d(np.arange(ngdl), c)  # d = [4 5 6 7 8 9] - 1

# %% extraigo las submatrices y especifico las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]; fd = f[c]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]; fc = f[d]

# desplazamientos para los gdls c = [1 2 3 10 11 12]
ac = np.array([0, 0, 0, 0, 0, 0]) # desplazamientos conocidos

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # calculo desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad                   # calculo fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl)  # separo la memoria
a[c] = ac;          a[d] = ad           # desplazamientos 
q[c] = qd         # q[d] = qc = 0       # fuerzas nodales de equilibrio

print('Desplazamientos de los nodos en coord. globales = ')
print(a)

# %% imprimo las fuerzas internas en cada barra referidas a las coordenadas
#    globales
for e in range(nb): # para cada barra
    print(f'\n\nFuerzas internas para barra {e+1} en coordenadas globales = ')
    qe_coord_glob = Ke[e]@a[idx[e,:]] - fe[e]
    print(qe_coord_glob)
    
    print(f'\nFuerzas internas para barra {e+1} en coordenadas locales =')
    print(T[e]@qe_coord_glob)
