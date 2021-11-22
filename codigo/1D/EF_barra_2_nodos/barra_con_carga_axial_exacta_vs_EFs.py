# -*- coding: utf-8 -*-

# %% definición del problema
# Calcule los desplazamientos y las reacciones en el empotramiento
# de la barra mostrada
#
# | b (carga distribuída de magnitud b)
# |->->->->->->->->->->->->->->->->
# |====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
# 1    2    3    4          nno-1  nno
# |<----longitud L de la barra---->|   el area transversal de la barra es A

import numpy as np
import matplotlib.pyplot as plt
# %matplotlib --list
# %matplotlib auto
# %matplotlib inline

# %% defino las variables
nef  = 3                       # número de elementos finitos (EF)
nno  = nef+1                   # número de nodos
ngdl = nno                     # número de grados de libertad
E    = 200e9     # Pa          # módulo de elasticidad de la barra
A    = (0.01)**2 # m^2         # área transversal de la barra
L    = 2         # m           # longitud de la barra
b    = 1000      # N/m         # fuerza axial aplicada sobre cada EF
P    = 250       # N           # carga nodal al final de la barra

xnod = np.linspace(0, L, nno)  # posición de los nodos
Le   = np.diff(xnod)           # longitud de cada EF (= np.tile(L/nef, (nef, 1))
k    = E*A/Le                  # rigidez de cada EF

# definición de EFs con respecto a nodos
rango = lambda a,b : np.arange(a,b+1) # np.arange con el punto final
LaG = np.column_stack((rango(1, nno-1), rango(2, nno))) - 1

# %% Relación de cargas puntuales
f = np.zeros(ngdl) # vector de fuerzas nodales equivalentes global
f[nno - 1] = P     # relaciono la carga puntual en el nodo "nno"

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K = np.zeros((ngdl,ngdl))   # matriz de rigidez global
for e in range(nef):        # ciclo sobre todos los elementos finitos
    idx = LaG[e,:]
    K[np.ix_(idx,idx)] += k[e]*np.array([[1., -1.], [-1., 1.]])
    f[idx]             += (b*Le[e]/2)*np.array([1., 1.])

# %% grados de libertad del desplazamiento conocidos y desconocidos
c = np.array([1]) - 1
d = np.setdiff1d(np.arange(ngdl), c)

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |
#|    | = |         ||    | - |    |     Recuerde que qc=0 (siempre)
#| qc |   | Kdc Kdd || ad |   | fc |

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[np.ix_(c,c)];  Kcd = K[np.ix_(c,d)]; fd = f[c]
Kdc = K[np.ix_(d,c)];  Kdd = K[np.ix_(d,d)]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = np.array([0])               # desplazamientos conocidos (en el gdl 1)

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # calculo desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # calculo fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # separo la memoria
a[c] = ac;          a[d] = ad          # desplazamientos
q[c] = qd         # q[d] = qc = 0      # fuerzas nodales de equilibrio

# %% calculo las cargas axiales en cada elemento finito
faxial = np.zeros(nef)
for e in range(nef):                   # ciclo sobre todas los elementos finitos
    Be = np.array([[-1/Le[e], 1/Le[e]]])
    ae = a[LaG[e,:]]
    faxial[e] = (E*A)*Be@ae            # = D*B(e)*a(e)

# %% imprimo los resultados
print('Desplazamientos (m) = \n',                             a[:,np.newaxis])
print('Fuerzas nodales equivalentes(N) = \n',                 f[:,np.newaxis])
print('Fuerzas nodales de equilibrio (N) = \n',               q[:,np.newaxis])
print('Cargas axiales en cada elemento finito (N) = \n', faxial[:,np.newaxis])

# %% Grafico la solución analítica y la solución por el MEF
# 1) grafico los desplazamientos de la barra
u = lambda x : (-b*x**2/2 + (P + b*L)*x)/(E*A) # solución analítica para el despl.

plt.figure                        # cree un nuevo lienzo
plt.subplot(2,1,1)                # grafique en la parte superior (1) del lienzo
xx = np.linspace(0, L, 100)       # 100 puntos equidistantes entre 0 y L
plt.plot(xx, u(xx), 'r', label='solución exacta de $u(x)$')
plt.plot(xnod, a, 'b.-', label='solución por el MEF')
plt.title('Comparación de la solución analítica con el MEF para el desplazamiento')
plt.xlabel('Eje X (m)')           # titulo del eje X
plt.ylabel('Desplazamiento (m)')  # titulo del eje Y
plt.legend(loc='lower right')

# 2) grafico la carga axial de la barra
faxial_exacta = lambda x : P + b*(L-x) # solución analítica para la carga axial

plt.subplot(2,1,2)                # grafique en la parte inferior (2) del lienzo
plt.plot(xx, faxial_exacta(xx), 'r', label='solución exacta de $f_{axial}(x)$')
for e in np.arange(nef):          # ciclo sobre todas los elementos finitos
    plt.plot([xnod[e], xnod[e+1]], [faxial[e], faxial[e]], 'b.-',
                                  label='solución por el MEF' if e == 0 else "")

plt.title('Comparación de la solución analítica con el MEF para la carga axial')
plt.xlabel('Eje X (m)')           # título del eje X
plt.ylabel('Carga axial (N)')     # título del eje Y
plt.legend(loc='upper right')
plt.show()

#bye, bye!!!
