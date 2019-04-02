#%% Ejemplo 11.3 Uribe Escamilla 

import numpy as np

# Unidades en tonedadas y cm

#%% defino las variables
X = 0
Y = 1
ang   = np.degrees(np.arctan2(300,400)) # angulo especificado en grados

#barra              1    2     3    4    5
theta = np.array([ang,   0, -ang,   0, -90 ])*np.pi/180 # angulo de inclinacion 
long  = np.array([500, 400,  500, 400, 300 ])           # longitud barra
area  = np.array([100,  40,  150,  40,  30 ])           # area barra

# LaG: local a global: matriz que relaciona nodos locales y globales
LaG = np.array([[1, 3],   # (se lee la barra x va del nodo i al nodo j)
                [1, 4],   # fila = barra
                [3, 2],   # col1 = nodo global asociado a nodo local 1
                [4, 2],   # col2 = nodo global asociado a nodo local 2
                [3, 4]]) - 1  

# gdl: grados de libertad
gdl = np.array([[1, 2],  # fila = nodo
                [3, 4],  # col1 = gdl en direccion x
                [5, 6],  # col2 = gdl en direccion y
                [7, 8]]) - 1  
      
# propiedades del material
E = 2040         # ton/cm^2
k = E*area/long  # rigidez de cada barra

#%% ensamblo la matriz de rigidez global
K = np.zeros((8,8))
T = 5*[None]        # MATLAB = cell(5,1) -> separo memoria
for e in range(5):  # para cada barra
    idx = np.r_[gdl[LaG[e,X],:], gdl[LaG[e,Y],:]]   # saco los 4 gdls de la barra
    c = np.cos(theta[e]); s = np.sin(theta[e]) # sin y cos de la inclinacion
    T[e] = np.array([[ c,  s,  0,  0],         # matriz de transformacion de coordenadas
                     [-s,  c,  0,  0],         # para la barra e
                     [ 0,  0,  c,  s], 
                     [ 0,  0, -s,  c]])
    Kloc = k[e]*np.array([[ 1,  0, -1,  0],    # matriz de rigidez local expresada en el 
                          [ 0,  0,  0,  0],    # sistema de coordenadas locales para la                  
                          [-1,  0,  1,  0],    # barra e
                          [ 0,  0,  0,  0]])
    K[np.ix_(idx,idx)] += T[e].T@Kloc@T[e] # sumo a K global

# %% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = np.array([1, 2, 4])-1;    d = np.array([3, 5, 6, 7, 8])-1

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]
# desplazamientos para los gdls c = [1 2 4]
ac = np.array([[0], [0], [0]])

# fuerzas en los gdls d = [3 5 6 7 8]
fc = np.array([[0],[5*np.cos(ang*np.pi/180)], [5*np.sin(ang*np.pi/180)], [0], [-20]]) # ton

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac)
qd = Kcc@ac + Kcd@ad

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros((8,1)); q = np.zeros((8,1))  # separo la memoria
a[c] = ac;     q[c] = qd
a[d] = ad    # q[d] = qc = 0

# %% calculo las cargas axiales (N) en cada barra
N = np.zeros(5)
for e in range(5): # para cada barra
    idx = np.r_[gdl[LaG[e,X],:], gdl[LaG[e,Y],:]] # saco los 4 gdls de la barra e
    N[e] = k[e]*np.array([-1, 0, 1, 0])@T[e]@a[idx]

#%% imprimo los resultados
print('a = \n', a, '\n')
print('q = \n', q, '\n')
print('N = \n', N, '\n')
