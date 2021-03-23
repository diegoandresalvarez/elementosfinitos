# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from misfunciones import dibujar_deformada, calc_Te, calc_feloc, calc_Keloc

#%% constantes que ayudarán a la lectura del código
NL1, NL2 = 0, 1
X, Y, TH = 0, 1, 2
g        = -9.81    # [m/s²] aceleración de la gravedad

# %% seleccione la malla a emplear:
nombre_archivo = 'fink'
#nombre_archivo = 'fink_portico'
#nombre_archivo = 'fink_cercha'
#nombre_archivo = 'torre_electrica'
#nombre_archivo = 'cercha_UribeEscamilla_11_3'
#nombre_archivo = 'cercha_UribeEscamilla_11_3_apoyo_inclinado'
#nombre_archivo = 'portico_UribeEscamilla_11_23'
df = pd.read_excel(f"{nombre_archivo}.xlsx", sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=0 o Y=1
xnod = df['xnod'][['x', 'y']].to_numpy()
nno  = xnod.shape[0]    # número de nodos (número de filas de la matriz xnod)

# %% definición de los grados de libertad
ngdl = 3*nno            # número de grados de libertad por nodo = [X, Y, TH]
# fila = nodo
# col1 = gdl en dirección x
# col2 = gdl en dirección y
# col3 = gdl en dirección angular antihoraria
gdl  = np.reshape(np.arange(ngdl), (nno, 3)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG = df['LaG_mat'][['NL1', 'NL2']].to_numpy() - 1
nbar = LaG.shape[0]      # número de EFs (número de filas de la matriz LaG)

# %% definición de los materiales
mat  = df['LaG_mat']['material'].to_numpy() - 1
E    = df['prop_mat']['E'].to_numpy()       # [Pa]     módulo de elasticidad
A    = df['prop_mat']['A'].to_numpy()       # [m²]     área transversal
I    = df['prop_mat']['I'].to_numpy()       # [m⁴]     inercia_y
rho  = df['prop_mat']['rho'].to_numpy()     # [kg/m³]  densidad
nmat = E.shape[0]                           # número de materiales

# %% definición del tipo de la barra
# 'EE'==pórtico, 'RR'==cercha, 'RE'/'ER' rotula-empotrado
tipo = df['LaG_mat']['tipo']

# %% relación de cargas puntuales
cp  = df['carga_punt']
ncp = cp.shape[0]       # número de cargas puntuales
f   = np.zeros(ngdl)    # vector de fuerzas nodales equivalentes global
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

#%% Se dibuja la estructura junto con su numeración
plt.figure(1)
for e in range(nbar):
   plt.plot(xnod[LaG[e,:],X], xnod[LaG[e,:],Y], 'b-')
   
   # Calculo la posición del centro de gravedad de la barra
   cgx = (xnod[LaG[e,NL1],X] + xnod[LaG[e,NL2],X])/2
   cgy = (xnod[LaG[e,NL1],Y] + xnod[LaG[e,NL2],Y])/2
   plt.text(cgx, cgy, str(e+1), color='red')

plt.plot(xnod[:,X], xnod[:,Y], 'ro')
for n in range(nno):
    plt.text(xnod[n,X], xnod[n,Y], str(n+1))
    
plt.axis('equal')
plt.grid(b=True, which='both', color='0.65',linestyle='-')
plt.title('Numeración de los nodos y barras de la estructura')
plt.show()

#%% fuerzas distribuidas aplicadas sobre las barras en coordenadas locales
b1  = df['carga_distr']['b1'].astype('float64').to_numpy()
b2  = df['carga_distr']['b2'].astype('float64').to_numpy()
q1  = df['carga_distr']['q1'].astype('float64').to_numpy()
q2  = df['carga_distr']['q2'].astype('float64').to_numpy()

# %% se ensambla en el vector de fuerzas nodales equivalentes global
K   = np.zeros((ngdl,ngdl))  # separo memoria
fe  = nbar*[None]
T   = nbar*[None]
L   = nbar*[None]
Ke  = nbar*[None]
idx = np.zeros((nbar,6), dtype=int)
for e in range(nbar):
   # saco los 6 gdls de la barra e
   idx[e] = np.r_[gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:]]
   
   # longitud de la barra
   x1, y1 = xnod[LaG[e,NL1], :]
   x2, y2 = xnod[LaG[e,NL2], :]
   L[e]  = np.hypot(x2-x1, y2-y1)
   
   # matriz de transformación Te
   T[e] = calc_Te(x1,y1, x2,y2)
   
   # vector fuerzas nodales equiv. y matriz de rigidez en coordenadas locales
   # incluye fuerzas por peso propio
   wx = rho[mat[e]]*A[mat[e]]*g * (y2-y1)/L[e]
   b1[e] += wx
   b2[e] += wx
   wy = rho[mat[e]]*A[mat[e]]*g * (x2-x1)/L[e]
   q1[e] += wy
   q2[e] += wy   
   
   feloc = calc_feloc(tipo[e], L[e], b1[e], b2[e], q1[e], q2[e])
   Keloc = calc_Keloc(tipo[e], L[e], A[mat[e]], E[mat[e]], I[mat[e]])
   
   # se convierten a coordenadas globales
   fe[e] = T[e].T @ feloc
   Ke[e] = T[e].T @ Keloc @ T[e]
   
   # se ensamblan
   f[idx[e]]                += fe[e] # ensambla fe{e} en f global
   K[np.ix_(idx[e],idx[e])] += Ke[e] # ensambla Ke{e} en K global

# %% restricciones y los grados de libertad del desplazamiento conocidos (c)
restric = df['restric']
nres    = restric.shape[0]
c       = np.empty(nres, dtype=int)
for i in range(nres):
   c[i] = gdl[restric['nodo'][i]-1, restric['dirección'][i]-1]

# desplazamientos conocidos
ac = restric['desplazamiento'].to_numpy()

# grados de libertad del desplazamiento desconocidos
d = np.setdiff1d(np.arange(ngdl), c)

# %% introduciendo los soportes inclinados
# rotaciones de los apoyos
ang = np.radians(restric['rotación'].to_numpy())

T_apoyo = np.identity(ngdl)   # FALTA HACERLA SPARSE
for i in range(nres):
    if not np.isnan(ang[i]):
        idx_res = gdl[restric['nodo'][i]-1, [X, Y]]
        T_apoyo[np.ix_(idx_res, idx_res)] = [[ np.cos(ang[i]), np.sin(ang[i])],
                                             [-np.sin(ang[i]), np.cos(ang[i])]]
    
# convierto a sistema de coordenadas con soportes inclinados 
# (OJO se sobreescriben K y f)
K = T_apoyo @ K @ T_apoyo.T
f = T_apoyo @ f

#%%
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    Recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[np.ix_(c,c)];  Kcd = K[np.ix_(c,d)]; fd = f[c]
Kdc = K[np.ix_(d,c)];  Kdd = K[np.ix_(d,d)]; fc = f[d]

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac)   # desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd                # fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); a[c] = ac; a[d] = ad # desplazamientos
q = np.zeros(ngdl); q[c] = qd            # fuerzas nodales equivalentes

# %% retorno las fuerzas y los desplazamientos en el sistema de coordenadas
#    donde los grados de libertad son paralelos a los ejes
# (OJO se sobreescriben q y a)
q = T_apoyo.T @ q
a = T_apoyo.T @ a

#%% calculo las fuerzas internas en cada barra qe
qe_loca = np.empty((nbar,6))
qe_glob = np.empty((nbar,6))
for e in range(nbar): # para cada barra
   qe_glob[e,:] = Ke[e]@a[idx[e]] - fe[e] # qe en coordenadas globales
   qe_loca[e,:] = T[e]@qe_glob[e]         # qe en coordenadas locales

# haga los números super pequeños iguales a cero, para mejorar la claridad   
qe_loca[np.abs(qe_loca) < 1e-12] = 0

# %% se configuran las unidades y las escalas para la graficación
# se configuran las unidades para los reportes
config = df['config'].set_index('variable')
U_FUER     = config.loc['U_FUER']['valor']
U_LONG     = config.loc['U_LONG']['valor']

# se configuran las escalas para la graficación
esc_def    = config.loc['esc_def']['valor']    # deformada
esc_faxial = config.loc['esc_faxial']['valor'] # diagrama de fuerzas axiales
esc_V      = config.loc['esc_V']['valor']      # diagrama de fuerzas cortantes
esc_M      = config.loc['esc_M']['valor']      # diagrama de momentos flectores

#%% Dibujar la estructura y su deformada
# se crean tablas para reportar los resultados nodales de: desplazamientos (a),
# fuerzas nodales equivalentes (f) y fuerzas nodales de equilibrio (q)
vect_mov = a.reshape((nno,3)) # vector de movimientos
q_global = q.reshape((nno,3))

xdef = xnod + esc_def*vect_mov[:,[X, Y]]

plt.figure(2)  
plt.title(f'Deformada (escalada {esc_def} veces)')
plt.xlabel(f'x [{U_LONG}]')
plt.ylabel(f'y [{U_LONG}]')
plt.axis('equal')

plt.figure(3)  
plt.title(f'Fuerza axial [{U_FUER}]')
plt.xlabel(f'x [{U_LONG}]')
plt.ylabel(f'y [{U_LONG}]')
plt.axis('equal')

plt.figure(4)  
plt.title(f'Fuerza cortante [{U_FUER}]')
plt.xlabel(f'x [{U_LONG}]')
plt.ylabel(f'y [{U_LONG}]')
plt.axis('equal')

plt.figure(5)  
plt.title(f'Momento flector [{U_FUER} {U_LONG}]')
plt.xlabel(f'x [{U_LONG}]')
plt.ylabel(f'y [{U_LONG}]')
plt.axis('equal')

for e in range(nbar):
   x1, y1  = xnod[LaG[e,NL1], :]
   x2, y2  = xnod[LaG[e,NL2], :]
   ae_loca = T[e] @ a[idx[e]]

   dibujar_deformada(tipo[e],
      A[mat[e]], E[mat[e]], I[mat[e]],
      x1,y1, x2,y2, b1[e],b2[e], q1[e],q2[e], 
      qe_loca[e], ae_loca,
      esc_def, esc_faxial, esc_V, esc_M)

plt.show()

# %% Reporte de los resultados:
for i in range(nres):
    if np.isnan(ang[i]):
        vect_mov[restric['nodo'][i]-1, TH] = np.NaN  
        q_global[restric['nodo'][i]-1, TH] = np.NaN

tabla_aq = pd.DataFrame(
    data   = np.c_[vect_mov, q_global],
    index  = np.arange(1, nno+1),
    columns= [f'u [{U_LONG}]', f'v [{U_LONG}]', 'theta [rad]', 
              f'qx [{U_FUER}]', f'qy [{U_FUER}]', f'qm [{U_FUER} {U_LONG}]'])
tabla_aq.index.name = '# nodo'

# fuerzas internas en coordenadas globales
tabla_qeglob = pd.DataFrame(
      data = qe_glob,
   index   = np.arange(1, nbar+1),
   columns = [f'qx1 [{U_FUER}]', f'qy1 [{U_FUER}]', f'qm1 [{U_FUER} {U_LONG}]', 
              f'qx2 [{U_FUER}]', f'qy2 [{U_FUER}]', f'qm2 [{U_FUER} {U_LONG}]'])
tabla_qeglob.index.name = '# bar'

# fuerzas internas en coordenadas locales
tabla_qeloca = pd.DataFrame(
     data = qe_loca,
  index   = np.arange(1, nbar+1),
  columns = [f'qxloc1 [{U_FUER}]', f'qyloc1 [{U_FUER}]', f'qm1 [{U_FUER} {U_LONG}]', 
             f'qxloc2 [{U_FUER}]', f'qyloc2 [{U_FUER}]', f'qm2 [{U_FUER} {U_LONG}]'])
tabla_qeloca.index.name = '# bar'

# se crea un archivo de MS EXCEL
archivo_resultados = f"resultados_{nombre_archivo}.xlsx"
writer = pd.ExcelWriter(archivo_resultados, engine = 'xlsxwriter')

# cada tabla hecha previamente es guardada en una hoja del archivo de Excel
tabla_aq.to_excel    (writer, sheet_name = 'aq', na_rep='NaN')
tabla_qeglob.to_excel(writer, sheet_name = 'qeglob')
tabla_qeloca.to_excel(writer, sheet_name = 'qeloca')
writer.save()

print(f'Cálculo finalizado. En "{archivo_resultados}" se guardaron los resultados.')

'''
#%% imprimo los resultados
print(' ')
print('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)')
print('~'*80)
qq = q.reshape((nno,3))
for i in range(nno):
   if not np.allclose(qq[i,:], np.array([0, 0, 0])):
      print('Nodo %3d qx = %12.4g ton, qy = %12.4g ton, mom = %12.4g ton*m' % 
         (i+1, qq[i,X], qq[i,Y], qq[i,TH]))
'''      
#%% bye, bye!
