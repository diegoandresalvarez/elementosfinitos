# -*- coding: utf-8 -*-

#%% Ejemplo 11.3 Uribe Escamilla

# Unidades en tonedadas y cm

import numpy as np

# define sin() y cos(), pero con un argumento en grados (como en MATLAB)
np.cosd = lambda x : np.cos(np.deg2rad(x))
np.sind = lambda x : np.sin(np.deg2rad(x))

# se definen algunas constantes que hacen el código más legible
NL1, NL2 = 0, 1

#%% defino las variables
ang   = np.degrees(np.arctan2(300,400)) # ángulo especificado en grados

#barra              1    2     3    4    5
theta = np.array([ang,   0, -ang,   0, -90 ]) # angulo de inclinacion
L     = np.array([500, 400,  500, 400, 300 ]) # Litud barra
A     = np.array([100,  40,  150,  40,  30 ]) # A barra

# LaG: local a global: matriz que relaciona nodos locales y globales
LaG = np.array([[1, 3],   # (se lee la barra x va del nodo i al nodo j)
                [1, 4],   # fila = barra
                [3, 2],   # col1 = nodo global asociado a nodo local 1
                [4, 2],   # col2 = nodo global asociado a nodo local 2
                [3, 4]]) - 1

# gdl: grados de libertad
gdl = np.array([[1, 2],  # fila = nodo
                [3, 4],  # col1 = gdl en dirección x
                [5, 6],  # col2 = gdl en dirección y
                [7, 8]]) - 1

# propiedades del material
E = 2040   # ton/cm^2
k = E*A/L  # rigidez de cada barra

#%% ensamblo la matriz de rigidez global
K   = np.zeros((8,8))
T   = 5*[None]        # MATLAB = cell(5,1) -> separo memoria
idx = 5*[None]
for e in range(5):  # para cada barra
    # saco los 4 gdls de la barra
    idx[e] = np.r_[gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:]]

    # matriz de transformacion de coordenadas para la barra e
    c = np.cosd(theta[e]); s = np.sind(theta[e]) # sin y cos de la inclinación
    T[e] = np.array([[ c,  s,  0,  0],           
                     [-s,  c,  0,  0],           
                     [ 0,  0,  c,  s],
                     [ 0,  0, -s,  c]])

    # matriz de rigidez local expresada en el sistema de coordenadas locales 
    # para la barra e
    Kloc = np.array([[ k[e],  0, -k[e],  0],    
                     [ 0,     0,  0,     0],    
                     [-k[e],  0,  k[e],  0],    
                     [ 0,     0,  0,     0]])

    # sumo a K global
    K[np.ix_(idx[e],idx[e])] += T[e].T@Kloc@T[e]

# %% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = np.array([1, 2, 4])-1;    d = np.array([3, 5, 6, 7, 8])-1

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[np.ix_(c,c)];  Kcd = K[np.ix_(c,d)]
Kdc = K[np.ix_(d,c)];  Kdd = K[np.ix_(d,d)]

# desplazamientos para los gdls c = [1 2 4]
ac = np.array([0, 0, 0])

# fuerzas en los gdls d = [3 5 6 7 8]
fc = np.array([0, 5*np.cosd(ang), 5*np.sind(ang), 0, -20]) # ton

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac)
qd = Kcc@ac + Kcd@ad

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(8); q = np.zeros(8)  # separo la memoria
a[c] = ac;       q[c] = qd
a[d] = ad      # q[d] = qc = 0

# %% calculo las fuerzas axiales (fax) en cada barra
fax = np.zeros(5)
for e in range(5): # para cada barra
    fax[e] = np.array([-k[e], 0, k[e], 0])@T[e]@a[idx[e]]

#%% imprimo los resultados
print('a = \n', a,   '\n')
print('q = \n', q,   '\n')
print('N = \n', fax, '\n')
