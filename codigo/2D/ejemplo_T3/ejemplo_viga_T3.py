# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import func_EF_T3
from func_EF_T3 import t2ft_T3, plot_esf_def

# %% CÁLCULO DE UNA VIGA CON ELEMENTOS FINITOS TRIANGULARES PARA TENSION PLANA
# Definición del problema
# Calcule los desplazamientos y las reacciones en los empotramiento, las
# deformaciones y los esfuerzos de la estructura en TENSION PLANA mostrada
# en la figura adjunta

# %% constantes que ayudaran en la lectura del código
X, Y          = 0, 1
NL1, NL2, NL3 = 0, 1, 2

# %% defino las variables/constantes del sólido
Ee   = 200e9 # [Pa] módulo de elasticidad del sólido
nue  = 0.30  # coeficiente de Poisson
te   = 0.10  # [m] espesor del sólido
rhoe = 7850. # [kg/m³] densidad
g    = 9.81  # [m/s²] aceleración de la gravedad

# %% Seleccione la malla a emplear
# 1) Malla del ejemplo de la clase
#df   = pd.read_excel('malla_ejemplo.xlsx', sheet_name=None)

# 2) Malla refinada (malla elaborada por David Felipe Cano Perdomo)
df = pd.read_excel('malla_refinada.xlsx', sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=1 o Y=2
xnod = df['xnod'][['x','y']].to_numpy()
nno  = xnod.shape[0]    # número de nodos (número de filas de xnod)

# %% definición de los grados de libertad
ngdl = 2*nno            # número de grados de libertad (dos por nodo)
gdl  = np.reshape(np.arange(ngdl), (nno,2)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG = df['LaG'][['NL1','NL2','NL3']].to_numpy() - 1
nef = LaG.shape[0]      # número de EFs (número de filas de LaG)

# %% Relación de cargas puntuales
cp  = df['carga_punt']
ncp = cp.shape[0]       # número de cargas puntuales
f   = np.zeros(ngdl)    # vector de fuerzas nodales equivalentes global
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

# %% Se dibuja la malla de elementos finitos
plt.figure()
cg = np.zeros((nef,2))  # almacena el centro de gravedad
for e in range(nef):
   idx_NL = [NL1, NL2, NL3, NL1]
   plt.plot(xnod[LaG[e, idx_NL], X], xnod[LaG[e, idx_NL], Y], 'b')

   # Calculo la posición del centro de gravedad del triángulo
   cg[e] = np.mean(xnod[LaG[e,:], :], axis=0)
   plt.text(cg[e,X], cg[e,Y], f'{e+1}', horizontalalignment='center',
                                        verticalalignment='center',   color='b')

plt.plot(xnod[:,X], xnod[:,Y], 'r*')
for i in range(nno):
   plt.text(xnod[i,X], xnod[i,Y], f'{i+1}', color='r')
plt.gca().set_aspect('equal', adjustable='box')   
plt.tight_layout()
plt.title('Malla de elementos finitos')
plt.show()

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
# equivalentes global
K   = np.zeros((ngdl,ngdl)) # matriz de rigidez global como RALA (sparse)                                                         SPARSE
B   = nef * [None]          # contenedor para las matrices de deformación
idx = nef * [None]          # indices asociados a los gdl del EF e

# matriz constitutiva del elemento para TENSION PLANA
De = np.array([[ Ee/(1-nue**2)    , Ee*nue/(1-nue**2),  0              ],
               [ Ee*nue/(1-nue**2), Ee/(1-nue**2)    ,  0              ],
               [ 0                , 0                ,  Ee/(2*(1+nue)) ]])

for e in range(nef):        # ciclo sobre todos los elementos finitos
   # Calculo de la matriz de rigidez del elemento e
   x1, y1 = xnod[LaG[e,NL1], :]
   x2, y2 = xnod[LaG[e,NL2], :]
   x3, y3 = xnod[LaG[e,NL3], :]

   Ae = 0.5*np.linalg.det(np.array([[ 1, x1, y1 ],      # área del EF e
                                    [ 1, x2, y2 ],
                                    [ 1, x3, y3 ]]))
   if Ae <= 0:
      raise Exception(
         f'La numeración local del EF {e+1} deben especificarse en sentido antihorario.\n')

   # Calculo de la matriz de deformaciones B.
   a1 = x2*y3 - x3*y2;        b1 = y2-y3;        c1 = x3-x2
   a2 = x3*y1 - x1*y3;        b2 = y3-y1;        c2 = x1-x3
   a3 = x1*y2 - x2*y1;        b3 = y1-y2;        c3 = x2-x1

   B[e] = (1/(2*Ae))*np.array([[ b1,  0,   b2,  0,   b3,  0 ],
                               [  0, c1,    0, c2,    0, c3 ],
                               [ c1, b1,   c2, b2,   c3, b3 ]])

   Ke = B[e].T@De@B[e]*te*Ae

   # Calculo del vector de fuerzas nodales equivalentes del elemento e
   # Fuerzas másicas (peso propio)
   fbe = np.array([0, -rhoe*g, 0, -rhoe*g, 0, -rhoe*g]) * Ae*te/3

   fe = fbe # vector de fuerzas nodales equivalentes

   # Ensamblo las contribuciones a las matrices globales
   idx[e] = np.r_[ gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:], gdl[LaG[e,NL3],:] ]
   K[np.ix_(idx[e],idx[e])] += Ke
   f[np.ix_(idx[e])]        += fe

# %% Muestro la configuración de la matriz K (K es rala)
plt.figure()
plt.spy(K)
plt.title('Los puntos representan los elementos diferentes de cero')
plt.show()

# %% Relación de las cargas superficiales (vector ft)
cd   = df['carga_distr']
nlcd = cd.shape[0]    # número de lados con carga distribuída
ft   = np.zeros(ngdl) # fuerzas nodales equivalentes de cargas superficiales                                             SPARSE
for i in range(nlcd):
   e     = cd['elemento'][i] - 1
   lado  = cd['lado'][i]
   carga = cd[['tix','tiy','tjx','tjy']].loc[i].to_numpy()
   fte = t2ft_T3(xnod[LaG[e,:],:], lado, carga, te)

   ft[np.ix_(idx[e])] += fte

# %% agrego al vector de fuerzas nodales equivalentes las fuerzas
# superficiales calculadas
f += ft

# %% restricciones y los grados de libertad del desplazamiento conocidos (c)
restric = df['restric']
nres = restric.shape[0]
c    = np.empty(nres, dtype=int)
for i in range(nres):
   c[i] = gdl[restric['nodo'][i]-1, restric['dirección'][i]-1]

# desplazamientos conocidos
ac = restric['desplazamiento'].to_numpy()

# grados de libertad del desplazamiento desconocidos
d = np.setdiff1d(range(ngdl), c)

# %% extraigo las submatrices y especifico las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |  # recuerde que qc=0 (siempre)
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |
Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]; fd = f[c]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]; fc = f[d]

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # separo la memoria
a[c] = ac;          a[d] = ad          # desplazamientos
q[c] = qd         # q[d] = qc = 0      # fuerzas nodales de equilibrio

# %% imprimo los resultados de los desplazamientos (a), las fuerzas nodales
# equivalentes (f) y nodales de equilibrio (q)
tabla_afq = pd.DataFrame(
    data=np.c_[a.reshape((nno,2)), f.reshape((nno,2)), q.reshape((nno,2))],
    index=np.arange(nno)+1,
    columns=['ux [m]', 'uy [m]', 'fx [N]', 'fy [N]', 'qx [N]', 'qy [N]'])
tabla_afq.index.name = '# nodo'

# %% Dibujo la malla de elementos finitos y las deformada de esta
delta  = np.reshape(a, (nno,2))
escala = 20000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

plt.figure()
for e in range(nef):
   nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]
   plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], 'r', 
                        label='Posición original'  if e == 0 else "", lw=0.5)
   plt.plot(xdef[nod_ef, X], xdef[nod_ef, Y], 'b', 
                        label='Posición deformada' if e == 0 else "")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.title(f'Deformada escalada {escala} veces')
plt.show()
plt.tight_layout()

# %% Se calcula para cada elemento las deformaciones y los esfuerzos
deform = np.zeros((3,nef))
esfuer = np.zeros((3,nef))
for e in range(nef):
   ae = a[idx[e]]               # desplazamientos de los gdl del elemento e
   deform[:,e] = B[e]@ae        # calculo las deformaciones
   esfuer[:,e] = De@deform[:,e] # calculo los esfuerzos

sx = esfuer[0,:];  sy = esfuer[1,:];  txy = esfuer[2,:]
ex = deform[0,:];  ey = deform[1,:];  gxy = deform[2,:]
ez = -(nue/Ee)*(sx+sy)          # deformaciones ez en tensión plana

# %% imprimo y grafico las deformaciones
tabla_exeyezgxy = pd.DataFrame(
   data=np.c_[ex, ey, ez, gxy],
   index=np.arange(nef)+1,
   columns=['ex', 'ey', 'ez', 'gxy [rad]'])
tabla_exeyezgxy.index.name = '# EF'

func_EF_T3.compartir_variables(xnod, LaG, cg)
plot_esf_def(ex,  r'$\epsilon_x$')
plot_esf_def(ey,  r'$\epsilon_y$')
plot_esf_def(ez,  r'$\epsilon_z$')
plot_esf_def(gxy, r'$\gamma_{xy}$')

# %% imprimo y grafico los esfuerzos
tabla_sxsytxy = pd.DataFrame(
   data=np.c_[sx, sy, txy],
   index=np.arange(nef)+1,
   columns=['sx [Pa]', 'sy [Pa]', 'txy [Pa]'])
tabla_sxsytxy.index.name = '# EF'

plot_esf_def(sx,  r'$\sigma_x$')
plot_esf_def(sy,  r'$\sigma_y$')
plot_esf_def(txy, r'$\tau_{xy}$')

# %% Se calculan y grafican para cada elemento los esfuerzos principales y
# sus direcciones
# NOTA: esto solo es válido para el caso de TENSION PLANA.
# En caso de DEFORMACIÓN PLANA se deben calcular los valores y vectores
# propios de la matriz de tensiones de Cauchy:
# [dirppales{e}, esfppales{e}] = eig([sx  txy 0    # matriz de esfuerzos
#                                     txy sy  0    # de Cauchy
#                                     0   0   0])

s1   = (sx+sy)/2 + np.sqrt(((sx-sy)/2)**2 + txy**2) # esfuerzo normal máximo
s2   = (sx+sy)/2 - np.sqrt(((sx-sy)/2)**2 + txy**2) # esfuerzo normal mínimo
tmax = (s1-s2)/2                                    # esfuerzo cortante máximo
ang  = 0.5*np.arctan2(2*txy, sx-sy)                 # ángulo asociado a s1

# %% Calculo de los esfuerzos de von Mises
s3 = np.zeros(nef)
sv = np.sqrt(((s1-s2)**2 + (s2-s3)**2 + (s1-s3)**2)/2)

# %% imprimo y grafico los esfuerzos s1, s2, tmax y sv
tabla_s1s2tmaxsv = pd.DataFrame(
   data=np.c_[s1, s2, tmax, sv, ang],
   index=np.arange(nef)+1,
   columns=['s1 [Pa]', 's2 [Pa]', 'tmax [Pa]', 'sv [Pa]', 'theta [rad]'])
tabla_s1s2tmaxsv.index.name = '# EF'

plot_esf_def(s1,   r'$\sigma_1$',   [ ang ])
plot_esf_def(s2,   r'$\sigma_2$',   [ ang+np.pi/2 ] )
plot_esf_def(tmax, r'$\tau_{max}$', [ ang-np.pi/4, ang+np.pi/4 ])
plot_esf_def(sv,   r'$\sigma_{VM}$')

# %%Se escriben los resultados a una hoja de MS EXCEL
nombre_archivo = 'resultados.xlsx'
writer = pd.ExcelWriter(nombre_archivo, engine='xlsxwriter')

# se escribe cada DataFrame a una hoja diferente
tabla_afq.to_excel(       writer, sheet_name='afq')
tabla_exeyezgxy.to_excel( writer, sheet_name='exeyezgxy')
tabla_sxsytxy.to_excel(   writer, sheet_name='sxsytxy')
tabla_s1s2tmaxsv.to_excel(writer, sheet_name='s1s2tmaxsv')

# Se cierra y graba el archivo de MS EXCEL
writer.save()
print(f'Los resultados se escribieron en {nombre_archivo}')

# %%bye, bye!