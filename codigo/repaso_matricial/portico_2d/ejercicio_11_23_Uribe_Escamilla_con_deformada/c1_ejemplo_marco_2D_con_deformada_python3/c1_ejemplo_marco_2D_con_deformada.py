# -*- coding: utf-8 -*-

#%% Unidades en kN y m

from math import atan2
import numpy as np
import matplotlib.pyplot as plt
from c1_dibujar_barra_deformada_portico import c1_dibujar_barra_deformada_portico

#%% constantes
LOC1 = 0; LOC2 = 1; X = 0; Y = 1; TH = 2

#%% Se define la estructura
xnod = np.array([[3, 4],   # coordenadas de cada nodo [x, y]
                 [7, 6],
                 [9, 0],
                 [0, 0]])

# LaG: local a global: matriz que relaciona nodos locales y globales
# fila = barra
# col1 = nodo global asociado a nodo local 1
# col2 = nodo global asociado a nodo local 2
# (se lee la barra x va del nodo i al nodo j)

#                 nod1  nod2 material
barra = np.array([[1,    2,    1],
                  [4,    1,    2],
                  [2,    3,    2]])-1

LaG = barra[:, [0, 1]]  # local a global
mat = barra[:, 2]       # material

#                  area       inercias_y       modulo de elasticidad
#                  A(m^2)     I(m^4)           E(ton/m^2)
props = np.array([[.30*.35,   .30*.35**3/12,     190e4],
                  [.30*.30,   .30*.30**3/12,     190e4]])

A = props[:,0];   I = props[:,1];   E = props[:,2]

nno  = xnod.shape[0] # numero de nodos (numero de filas de xnod)
nbar = LaG.shape[0]  # numero de EFs (numero de filas de LaG)
ngdl = 3*nno         # numero de grados de libertad (tres por nodo)

#%% gdl: grados de libertad
# fila = nodo
# col1 = gdl en dirección x
# col2 = gdl en dirección y
# col3 = gdl en dirección angular antihoraria
gdl = np.arange(ngdl).reshape(nno, 3)  # nodos vs grados de libertad

#%% cargas aplicadas (gdl carga)
cargas_aplica = np.array([[gdl[1-1,X], 1.5]])
dofs_cargados = cargas_aplica[:,0].astype(int)

f = np.zeros((ngdl, 1))
f[dofs_cargados] = cargas_aplica[:,1]

#%% Se dibuja la estructura junto con su numeración
plt.figure(1)
for e in range(nbar):
   plt.plot(xnod[LaG[e,:],X], xnod[LaG[e,:],Y], 'b-')
   
   # Calculo la posición del centro de gravedad de la barra
   cgx = (xnod[LaG[e,LOC1],X] + xnod[LaG[e,LOC2],X])/2
   cgy = (xnod[LaG[e,LOC1],Y] + xnod[LaG[e,LOC2],Y])/2
   plt.text(cgx, cgy, str(e+1), color='red')

plt.plot(xnod[:,X], xnod[:,Y], 'ro')
for n in range(nno):
    plt.text(xnod[n,X], xnod[n,Y], str(n+1))
    
plt.axis('equal')
plt.grid(b=True, which='both', color='0.65',linestyle='-')
plt.title('Numeración de la estructura')
plt.show()

#%% fuerzas distribuidas aplicadas sobre las barras en coordenadas locales
ang1 = atan2(2,4);

qxloc = [
   lambda x: -2.8*np.sin(ang1)*np.cos(ang1),
   lambda x: 0,
   lambda x: 0
]

qyloc = [
   lambda x: -2.8*np.cos(ang1)**2,
   lambda x: 0,
   lambda x: 0
]

#%% fuerzas nodales equivalentes para las diferentes barras
# (en este ejemplo las fuerzas nodales equivalentes estas siendo 
# especificadas con respecto al sistema de coordenadas globales)

fe = nbar*[None]
#                  fxi    fyi      mi     fxj      fyj     mj
#                  ton    ton      ton-m  ton      ton     ton-m
fe[0] = np.array([[0,     -5.60,  -3.733,  0,     -5.60,   +3.733 ]]).T # OJO con los signos
fe[1] = np.array([[0,      0   ,   0    ,  0,      0   ,   0      ]]).T # mirar pag 613
fe[2] = np.array([[0,      0   ,   0    ,  0,      0   ,   0      ]]).T

#%% ensamblo la matriz de rigidez global
K   = np.zeros((ngdl,ngdl))  # separo memoria
Ke  = nbar*[None]
T   = nbar*[None]
idx = np.zeros((nbar,6), dtype=int)
for e in range(nbar): # para cada barra
   # saco los 6 gdls de la barra e
   idx[e,:] = np.r_[gdl[LaG[e,LOC1],:], gdl[LaG[e,LOC2],:]]
   
   x1 = xnod[LaG[e,LOC1], X];  x2 = xnod[LaG[e,LOC2], X]
   y1 = xnod[LaG[e,LOC1], Y];  y2 = xnod[LaG[e,LOC2], Y]
   
   L = np.hypot(x2-x1, y2-y1)
   
   # matriz de transformación de coordenadas para la barra e
   c = (x2-x1)/L;   s = (y2-y1)/L;  # seno y coseno de la inclinación
   T[e] = np.array([[ c,  s,  0,  0,  0,  0],
                    [-s,  c,  0,  0,  0,  0],
                    [ 0,  0,  1,  0,  0,  0],
                    [ 0,  0,  0,  c,  s,  0],
                    [ 0,  0,  0, -s,  c,  0],
                    [ 0,  0,  0,  0,  0,  1]])
         
   # matriz de rigidez local expresada en el sistema de coordenadas locales
   # para la barra e
   AE = A[mat[e]]*E[mat[e]];       L2=L**2
   EI = E[mat[e]]*I[mat[e]];       L3=L**3
   Kloc = np.array([
        [ AE/L,   0      ,   0      ,  -AE/L,    0      ,   0      ],  
        [ 0   ,  12*EI/L3,   6*EI/L2,   0   ,  -12*EI/L3,   6*EI/L2],
        [ 0   ,   6*EI/L2,   4*EI/L ,   0   ,   -6*EI/L2,   2*EI/L ],
        [-AE/L,   0      ,   0      ,   AE/L,    0      ,   0      ],
        [ 0   , -12*EI/L3,  -6*EI/L2,   0   ,   12*EI/L3,  -6*EI/L2],
        [ 0   ,   6*EI/L2,   2*EI/L ,   0   ,   -6*EI/L2,   4*EI/L ]])

   # matriz de rigidez local en coordenadas globales
   Ke[e] = T[e].T @ Kloc @ T[e]
   K[np.ix_(idx[e,:],idx[e,:])] += Ke[e] # sumo Ke[e] a K global
   f[idx[e,:]]                  += fe[e] # sumo a f global

#%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
apoyos = np.array([[gdl[3-1,X],  0],
                   [gdl[3-1,Y],  0],
                   [gdl[3-1,TH], 0],
                   [gdl[4-1,X],  0],
                   [gdl[4-1,Y],  0],
                   [gdl[4-1,TH], 0]])

c = apoyos[:,0].astype(int)
d = np.setdiff1d(np.arange(ngdl), c)

#%%
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    Recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]; fd = f[c]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = apoyos[:,1][:,np.newaxis] # desplazamientos conocidos en contorno

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros((ngdl,1)); a[c] = ac; a[d] = ad # desplazamientos
q = np.zeros((ngdl,1)); q[c] = qd            # fuerzas nodales equivalentes

#%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
# globales
qe_coord_loc   = nbar*[None]
qe_coord_glob = nbar*[None]
for e in range(nbar): # para cada barra
   print('\n\n Fuerzas internas para barra %d en coord. globales =' % (e+1))
   qe_coord_glob[e] = Ke[e]@a[idx[e,:]] - fe[e]
   print(qe_coord_glob[e])
   
   print('\n\n Fuerzas internas para barra %d en coord. locales =' % (e+1))
   qe_coord_loc[e] = T[e]@qe_coord_glob[e]
   print(qe_coord_loc[e])

#%% imprimo los resultados
print('Desplazamientos nodales                                               ')
print('~'*80)
vect_mov = a.reshape((nno,3)) # vector de movimientos
for i in range(nno):
   print('Nodo %3d: u = %12.4g mm, v = %12.4g mm, theta = %12.4g rad' %
      (i+1, 1000*vect_mov[i,X], 1000*vect_mov[i,Y], vect_mov[i,TH]))

print(' ')
print('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)')
print('~'*80)
qq = q.reshape((nno,3))
for i in range(nno):
   if not np.allclose(qq[i,:], np.array([0, 0, 0])):
      print('Nodo %3d qx = %12.4g kN, qy = %12.4g kN, mom = %12.4g kN*m' % 
         (i+1, qq[i,X], qq[i,Y], qq[i,TH]))

#%% Dibujar la estructura y su deformada
esc_def    = 50            # escalamiento de la deformada
esc_faxial = 0.2           # escalamiento del diagrama de axiales
esc_V      = 0.3           # escalamiento del diagrama de cortantes
esc_M      = 0.3           # escalamiento del diagrama de momentos

xdef = xnod + esc_def*vect_mov[:,[X, Y]]

plt.figure(2)  
plt.title('Deformada')
plt.xlabel('x, m')
plt.ylabel('y, m')
plt.axis('equal')

plt.figure(3)  
plt.title('Fuerza axial [kN]')
plt.xlabel('x, m')
plt.ylabel('y, m')
plt.axis('equal')

plt.figure(4)  
plt.title('Fuerza cortante [kN]')
plt.xlabel('x, m')
plt.ylabel('y, m')
plt.axis('equal')

plt.figure(5)  
plt.title('Momento flector [kN-m]')
plt.xlabel('x, m')
plt.ylabel('y, m')
plt.axis('equal')

for e in range(nbar):
   x1 = xnod[LaG[e,LOC1], X];  x2 = xnod[LaG[e,LOC2], X]
   y1 = xnod[LaG[e,LOC1], Y];  y2 = xnod[LaG[e,LOC2], Y]

   c1_dibujar_barra_deformada_portico(
      A[mat[e]], E[mat[e]], I[mat[e]],
      x1,y1, x2,y2, qxloc[e], qyloc[e], 
      np.squeeze(qe_coord_loc[e]), np.squeeze(T[e] @ a[idx[e,:]]),
      esc_def, esc_faxial, esc_V, esc_M)

plt.show()

#%% bye, bye!
