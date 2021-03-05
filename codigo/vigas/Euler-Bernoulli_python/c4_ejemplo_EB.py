# -*- coding: utf-8 -*-

# Programa para el cálculo de vigas de Euler-Bernoulli.

# %%Importación de librerías

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path

# %%Defino las constantes y variables
Y = 0; TH = 1  # Y: vertical, TH: rotacional

#filename = 'viga_Uribe_Escamilla_ej_5_5'
filename = 'viga_con_resortes'
archivo_xlsx = pd.read_excel(path.join('..', 'ejemplos', f'{filename}.xlsx'),
                             sheet_name=None)

# %%Se lee la posición de los nodos
T       = archivo_xlsx['xnod']
idxNODO = np.array(T['nodo'])
xnod    = np.array(T['x'])                      # Posicion de los nodos
L       = np.diff(xnod)                         # Longitud de cada EF

nno  = len(xnod)                                # Número de nodos
nef  = nno - 1                                  # Número de element finitos(EF)
ngdl = 2*nno                                    # Número de grados de libertad
gdl  = np.array([range(ngdl)]).reshape(nno, 2)  # Grados de libertad

# %%Se leen la matriz de conectividad (LaG), el modulo de elasticidad, las 
# propiedades del material y las cargas
T     = archivo_xlsx['LaG_EI_q']
idxEF = np.array(T['EF']) - 1
LaG   = np.array(T[['NL1', 'NL2']]) -1  # Definición de EFs respecto a nodos
E     = np.array(T['E'])                # Módulo de elasticidad E del EF
I     = np.array(T['I'])                # Momento de inercia Iz del EF
G     = np.array(T['G'])                # Módulo de rigidez (para viga de Tim)
Aast  = np.array(T['Aast']);            # Área de cortante (para viga de Tim)
fz    = T[['q1e', 'q2e']]               # Relación de las cargas distribuidas
fz    = np.array(fz.fillna(0))          # Reemplazo de los NaN con ceros

# %%Relación de los apoyos
T       = archivo_xlsx['restric']
idxNODO = np.array(T['nodo'])            # Se quita 1 para el indexado de Py
dirdesp = np.array(T['direccion'])       # Se quita 1 para el indexado de Py
ac      = np.array(T['desplazamiento'])  # Desplazamientos conocidos

# %%Grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = len(idxNODO)
c = np.empty(n_apoyos, dtype=int)                 # GDL conocidos
for i in range(n_apoyos):
    c[i] = gdl[idxNODO[i] - 1, dirdesp[i] - 1]

d = np.setdiff1d(range(ngdl), c)  # GDL desconocidos

# %%Relación de cargas puntuales
T       = archivo_xlsx['carga_punt']
idxNODO = np.array(T['nodo'])       # Se quita 1 para el indexado de Py
dirfp   = np.array(T['direccion'])  # Se quita 1 para el indexado de Py
fp      = np.array(T['fuerza_puntual'])

# %%Se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales
# equivalentes global "f"
f_ini = np.zeros(ngdl)  # Vector de fuerzas nodales equivalentes global
for i in range(len(idxNODO)):
    f_ini[gdl[idxNODO[i] - 1, dirfp[i] - 1]] = fp[i]

# %%Relacion de los resortes
T       = archivo_xlsx['resortes']
idxNODO = np.array(T['nodo'])  # Se quita 1 para el indexado de Py
tipores = np.array(T['tipo'])  # Y=0 (vertical), TH=1 (rotacional)
kres    = np.array(T['k'])     # Constante del resorte

# %%Grados de libertad del desplazamiento conocidos y desconocidos
K_ini = np.zeros((ngdl, ngdl))  # Matriz de rigidez global
n_resortes = len(idxNODO)
for i in range(n_resortes):
    idx = gdl[idxNODO[i] - 1, tipores[i] - 1]
    K_ini[idx, idx] = kres[i]


# %%VIGA DE EULER-BERNOULLI:
# Con el programa "func_forma_euler_bernoulli.m" se calcularon:
#   Ke     = la matriz de rigidez de flexion del elemento e
#   fe     = el vector de fuerzas nodales equivalentes
#   Bb     = la matriz de deformaciones de flexion
#   N      = matriz de funciones de forma
#   dN_dxi = derivada de la matriz de funciones de forma con respecto a xi

# %%ensamblo la matriz de rigidez global y el vector de fuerzas nodales
# equivalentes global para la viga de Euler-Bernoulli
idx = nef * [None]  # grados de libertad del elemento e
K = K_ini
f = f_ini

for e in range(nef):     # Ciclo sobre todos los elementos finitos
    idx[e] = [gdl[LaG[e, 0], Y ],
              gdl[LaG[e, 0], TH],
              gdl[LaG[e, 1], Y ],
              gdl[LaG[e, 1], TH]]

    Le = L[e]

    # Matriz de rigidez de flexion del elemento e
    Ke = (E[e]*I[e]/Le**3) * np.array([[  12,    6*Le,   -12,    6*Le    ],
                                       [   6*Le, 4*Le**2, -6*Le, 2*Le**2 ],
                                       [ -12,   -6*Le,    12,   -6*Le    ],
                                       [   6*Le, 2*Le**2, -6*Le, 4*Le**2 ]])
    
    # Vector de fuerzas nodales equivalentes de una carga trapezoidal
    fe = np.array([ (   Le*(7*fz[e, 0] + 3*fz[e, 1]))/20,   # = Y1
                    (Le**2*(3*fz[e, 0] + 2*fz[e, 1]))/60,   # = M1
                    (   Le*(3*fz[e, 0] + 7*fz[e, 1]))/20,   # = Y2
                   -(Le**2*(2*fz[e, 0] + 3*fz[e, 1]))/60])  # = M2

    # Se ensambla la matriz de rigidez K y el vector de fuerzas nodales
    # equivalentes f
    K[np.ix_(idx[e], idx[e])] += Ke
    f[np.ix_(idx[e])]         += fe

# %%Se resuelve el sistema de ecuaciones
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd     |   | Kcc Kcd || ac |   | fd |
#|        | = |         ||    | - |    |
#| qc = 0 |   | Kdc Kdd || ad |   | fc |

# Extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[np.ix_(c, c)];  Kcd = K[np.ix_(c, d)]; fd = f[c]
Kdc = K[np.ix_(d, c)];  Kdd = K[np.ix_(d, d)]; fc = f[d]

ad = np.linalg.solve(Kdd, fc - Kdc@ac) # Desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # Fuerzas de equilibrio desconocidas

# Armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # Separo la memoria
a[c] = ac;          a[d] = ad          # Desplazamientos
q[c] = qd         # q[d] = qc = 0      # Fuerzas nodales equivalentes

# %%Cálculo de los momentos flectores y las fuerzas cortantes
# M = se calcula en las raices del polinomio de GL de orden 2
# V = se calcula en el centro del EF (raiz del polinomio de GL de orden 1)
# se reserva la memoria
xmom = np.zeros((2, nef))  # Posición donde se calcula
mom  = np.zeros((2, nef))  # Momento flector
cor  = np.zeros(nef)  # Fuerza cortante
xi =np.array([-np.sqrt(1/3), np.sqrt(1/3)]) # Raices del polinom de Legendre de grado 2

for e in range(nef):
    # Longitud del elemento finito e
    Le = L[e]

    # Matriz de deformaciones de flexión
    Bbe = np.array([ (6*xi)/Le**2,
                     (3*xi - 1)/Le,
                    -(6*xi)/Le**2,
                     (3*xi + 1)/Le]).T

    # Lugar donde se calcula el momento (centro del EF)
    xmom[:, e] = Le*xi/2 + (xnod[LaG[e, 0]] + xnod[LaG[e, 1]])/2

    # Vector de desplazamientos nodales del elemento ae
    ae = a[idx[e]]
    
    mom[:, e] = E[e]*I[e]*Bbe@ae                   # Momento flector
    dN3_dxi3  = np.array([ 3/2, (3*Le)/4, -3/2, (3*Le)/4 ]).T
    cor[e]    = E[e]*I[e]*(8/(Le**3))*dN3_dxi3@ae  # Fuerza cortante

# %%Se calculan los desplazamientos al interior de cada EF
nint = 10         # Número de puntos donde se interpolara dentro del EF
xi = np.linspace(-1,1,nint)  # Coordenadas naturales

xx = np.zeros((nef, nint))    # Interpol de posiciones (geometria) en el elemento
ww = np.zeros((nef, nint))    # Interpol desplazamientos en el elemento
tt = np.zeros((nef, nint))    # interpol angulo en el elemento
for e in range(nef):  # Ciclo sobre todas los elementos finitos
    # Longitud del elemento finito e
    Le = L[e]

    # Matriz de funciones de forma y su derivada
    N = np.array([xi**3/4 - (3*xi)/4 + 1/2,
                  -(Le*(- xi**3/4 + xi**2/4 + xi/4 - 1/4))/2,
                  - xi**3/4 + (3*xi)/4 + 1/2,
                  -(Le*(- xi**3/4 - xi**2/4 + xi/4 + 1/4))/2]).T

    dN_dxi = np.array([(3*xi**2)/4 - 3/4,
                        -(Le*(- (3*xi**2)/4 + xi/2 + 1/4))/2,
                        3/4 - (3*xi**2)/4,
                        (Le*((3*xi**2)/4 + xi/2 - 1/4))/2]).T

    # Vector de desplazamientos nodales del elemento a^{(e)}
    ae = a[idx[e]]

    # Interpola sobre la geometria (coord naturales a geometricas)
    xx[e] = Le*xi/2 + (xnod[LaG[e, 0]] + xnod[LaG[e, 1]])/2

    # Se calcula el desplazamiento al interior del elemento finito
    ww[e] = N@ae

    # Se calcula el angulo al interior del elemento finito
    tt[e] = np.arctan((dN_dxi*2/Le)@ae)

# %%Imprimo los resultados
print('Desplazamientos nodales                      ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
for n in range(nno):
    print('{3}Nodo{0:4d}:  w ={3} {3}{1:12.4} m, theta ={3} {3}{2:12.4} rad{3}'
          .format(n+1, a[2*n], a[2*n+1], ''))

print('\nFuerzas nodales de equilibrio '
      + '(solo se imprimen los diferentes de cero)')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
for n in range(nno):
    if not(q[2*n] == 0 and q[2*n+1] == 0):
        print('{3}Nodo{0:4d}:  Ry ={3} {3}{1:12.4} kN, Mz ={3} {3}{2:12.4} kN-m{3}'
              .format(n+1, q[2*n], q[2*n+1], ''))

# %%Gráfico de la solucion analitica y la solución por el MEF
# %%1) grafico los desplazamientos de la viga
fig, ax1 = plt.subplots(figsize=(8, 2.5))
for e in range(nef):  # Ciclo sobre todos los elementos finitos
    if e==0: label='Elementos finitos'
    else: label=None
    ax1.plot(xx[e], ww[e], 'b-', label=label, linewidth=1.2)

ax1.set_title('Solución con el MEF para el desplazamiento')
ax1.set_xlabel('Eje X (m)')           # Título del eje X
ax1.set_ylabel('Desplazamiento (m)')  # Título del eje Y
ax1.grid()                            # reticula
plt.xlim(xnod[1], xnod[-1])           # Rango en el eje X del gráfico

# %%2) Gráfico de los ángulos de giro
fig, ax2 = plt.subplots(figsize=(8, 2.5))
for e in range(nef):  # Ciclo sobre todos los elementos finitos
    if e==0: label='Elementos finitos' 
    else: label=None
    ax2.plot(xx[e], tt[e], 'b-', label=label, linewidth=1.2)

ax2.set_title('Solución con el MEF para el giro')
ax2.set_xlabel('Eje X (m)')   # Título del eje X
ax2.set_ylabel('Giro (rad)')  # Título del eje Y

ax2.grid()                    # reticula
plt.xlim(xnod[1], xnod[-1])   # Rango en el eje X del gráfico

# %%3) grafico los momentos
fig, ax3 = plt.subplots(figsize=(8, 2.5))
label='Elementos finitos'
ax3.plot(np.ravel(xmom, order='F'), np.ravel(mom, order='F'),
                'b-', label=label, linewidth=1.2)

ax3.set_title('Solucion con el MEF para el momento flector')

ax3.set_xlabel('Eje X (m)')               # Título del eje X
ax3.set_ylabel('Momento flector (kN-m)')  # Título del eje Y

ax3.grid()                                # reticula
plt.xlim(xnod[1], xnod[-1])               # Rango en el eje X del gráfico

# %%4) Gráfico de la fuerza cortante
fig, ax4 = plt.subplots(figsize=(8, 2.5))
for e in range(nef):  # Ciclo sobre todos los elementos finitos
    if e==0: label='Elementos finitos'
    else: label=None
    ax4.plot([xnod[LaG[e, 0]], xnod[LaG[e, 1]]],
                    [cor[e], cor[e]], 'b-', label=label, linewidth=1.2)

ax4.set_title('Solucion con el MEF para la fuerza cortante')

ax4.set_xlabel('Eje X (m)')             # Título del eje X
ax4.set_ylabel('Fuerza cortante (kN)')  # Título del eje Y

ax4.grid()                              # reticula
plt.xlim(xnod[1], xnod[-1])             # Rango en el eje X del gráfico

# %%Comparación con la solucion exacta (calculada con MAXIMA y el método de
# las funciones de discontinuidad

if filename == 'viga_con_resortes':
    import re
    def leer(fn):
        with open(
                path.join('..', 'ejemplos', 'results_viga_con_resortes_EB', fn),
                'r') as f:
            txt = f.read()
        nums = re.findall(r'\[([^][]+)\]', txt)
        return np.loadtxt(nums, delimiter=",")
   
    x = leer('x.txt')
    V = leer('Vx.txt')
    M = leer('Mx.txt')
    t = leer('tx.txt')
    w = leer('vxx.txt')

    ax1.plot(x, w, 'r.', label='Solución teórica', markersize=5)
    ax2.plot(x, t, 'r.', label='Solución teórica', markersize=5)
    ax3.plot(x, M, 'r.', label='Solución teórica', markersize=5)
    ax4.plot(x, V, 'r.', label='Solución teórica', markersize=5)

    ax1.legend(); ax2.legend(); ax3.legend(); ax4.legend()

# Fin del programa
