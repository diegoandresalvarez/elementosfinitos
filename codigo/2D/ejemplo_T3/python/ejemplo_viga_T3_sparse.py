# -*- coding: utf-8 -*-

# %%
'''
-------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para TENSION PLANA usando elementos
      triangulares de 3 nodos
-------------------------------------------------------------------------------

DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramiento, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func_EF_T3 import t2ft_T3, plot_esf_def, compartir_variables
from scipy import sparse 
from scipy.sparse.linalg import spsolve

# %% constantes que ayudarán en la lectura del código
X, Y          = 0, 1
NL1, NL2, NL3 = 0, 1, 2
g             = 9.81  # [m/s²]  aceleración de la gravedad

# %% Seleccione la malla a emplear
# 1) Malla del ejemplo de la clase
# nombre_archivo = 'malla1'            
# 2) Malla refinada (malla elaborada por David Felipe Cano Perdomo)
nombre_archivo = 'malla_refinada_v1' 
# 3) Malla extremadamente refinada cerca a las cargas puntuales y los apoyos
# nombre_archivo = 'malla_refinada_v2' 

df = pd.read_excel(f"{nombre_archivo}.xlsx", sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=0 o Y=1
xnod = df['xnod'][['x','y']].to_numpy()
nno  = xnod.shape[0]    # número de nodos (número de filas de la matriz xnod)

# %% definición de los grados de libertad
ngdl = 2*nno            # número de grados de libertad (dos por nodo)
gdl  = np.reshape(np.arange(ngdl), (nno,2)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG = df['LaG_mat'][['NL1','NL2','NL3']].to_numpy() - 1
nef = LaG.shape[0]      # número de EFs (número de filas de la matriz LaG)

# %% definición de los materiales
mat = df['LaG_mat']['material'].to_numpy() - 1
Ee   = df['prop_mat']['E'].to_numpy()       # [Pa]     módulo de elasticidad del sólido
nue  = df['prop_mat']['nu'].to_numpy()      # [-]      coeficiente de Poisson
rhoe = df['prop_mat']['rho'].to_numpy()     # [kg/m³]  densidad
te   = df['prop_mat']['espesor'].to_numpy() # [m]      espesor
nmat = Ee.shape[0]                          # número de materiales

# %% Relación de cargas puntuales
cp  = df['carga_punt']
ncp = cp.shape[0]       # número de cargas puntuales
f   = np.zeros(ngdl)    # vector de fuerzas nodales equivalentes global
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

# %% Se dibuja la malla de elementos finitos
cg = np.zeros((nef,2))  # almacena el centro de gravedad de los EF
plt.figure()
for e in range(nef):
   nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]
   plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], 'b')

   # se calcula la posición del centro de gravedad del triángulo
   cg[e] = np.mean(xnod[LaG[e,:], :], axis=0)

   # y se reporta el número del elemento actual
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
# K   = np.zeros((ngdl,ngdl))          # matriz de rigidez global
K   = sparse.coo_matrix((ngdl, ngdl)) # matriz de rigidez global como RALA (sparse)
B   = nef * [None]          # contenedor para las matrices de deformación
idx = nef * [None]          # indices asociados a los gdl del EF e

# matriz constitutiva del elemento para TENSION PLANA
De = nmat * [ None ]
for i in range(nmat):
    De[i] = np.array([[Ee[i]/(1 - nue[i]**2),        Ee[i]*nue[i]/(1 - nue[i]**2), 0                     ],
                      [Ee[i]*nue[i]/(1 - nue[i]**2), Ee[i]/(1 - nue[i]**2),        0                     ],
                      [0,                            0,                            Ee[i]/(2*(1 + nue[i]))]])

for e in range(nef):        # ciclo sobre todos los elementos finitos
   # Calculo de la matriz de rigidez del elemento e
   x1, y1 = xnod[LaG[e,NL1]]
   x2, y2 = xnod[LaG[e,NL2]]
   x3, y3 = xnod[LaG[e,NL3]]

   Ae = 0.5*np.linalg.det(np.array([[ 1, x1, y1 ],      # área del EF e
                                    [ 1, x2, y2 ],
                                    [ 1, x3, y3 ]]))
   if Ae <= 0:
      raise Exception(
         f'La numeración local del EF {e+1} debe especificarse en sentido antihorario.\n')

   # Calculo de la matriz de deformaciones B.
   a1 = x2*y3 - x3*y2;        b1 = y2-y3;        c1 = x3-x2
   a2 = x3*y1 - x1*y3;        b2 = y3-y1;        c2 = x1-x3
   a3 = x1*y2 - x2*y1;        b3 = y1-y2;        c3 = x2-x1

   B[e] = (1/(2*Ae))*np.array([[ b1,  0,   b2,  0,   b3,  0 ],
                               [  0, c1,    0, c2,    0, c3 ],
                               [ c1, b1,   c2, b2,   c3, b3 ]])

   Ke = te[mat[e]]*B[e].T@De[mat[e]]@B[e]*Ae

   # Calculo del vector de fuerzas nodales equivalentes del elemento e
   # Fuerzas másicas (peso propio)
   fbe = -rhoe[mat[e]]*g*Ae*te[mat[e]]*np.array([0., 1., 0., 1., 0., 1.])/3

   fe = fbe # vector de fuerzas nodales equivalentes

   # Ensamblo las contribuciones a las matrices globales
   idx[e] = np.r_[ gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:], gdl[LaG[e,NL3],:] ]
   # K[np.ix_(idx[e],idx[e])] += Ke
   IDX = np.array([(i,j) for i in idx[e] for j in idx[e]])
   K += sparse.coo_matrix((Ke.flat, (IDX[:,0],IDX[:,1])), shape=(ngdl,ngdl))
   f[np.ix_(idx[e])]        += fe

# %% Muestro la configuración de la matriz K (K es rala)
plt.figure()
plt.spy(K)
plt.title('Los puntos representan los elementos diferentes de cero')
plt.show()

# %% Relación de las cargas superficiales (vector ft)
cd   = df['carga_distr']
nlcd = cd.shape[0]    # número de lados con carga distribuída
ft   = np.zeros(ngdl) # fuerzas nodales equivalentes de cargas superficiales
for i in range(nlcd):
   e     = cd['elemento'][i] - 1
   lado  = cd['lado'][i]
   carga = cd[['tix','tiy','tjx','tjy']].loc[i].to_numpy()
   fte = t2ft_T3(xnod[LaG[e,:],:], lado, carga, te[mat[e]])

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
Kcc = K[np.ix_(c,c)];  Kcd = K[np.ix_(c,d)]; fd = f[c]
Kdc = K[np.ix_(d,c)];  Kdd = K[np.ix_(d,d)]; fc = f[d]

# %% resuelvo el sistema de ecuaciones
# ad = np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
ad = spsolve(Kdd, fc - Kdc@ac)
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
   ae = a[idx[e]]                       # desplazamientos de los gdl del elemento e
   deform[:,e] = B[e]@ae                # calculo las deformaciones
   esfuer[:,e] = De[mat[e]]@deform[:,e] # calculo los esfuerzos

sx = esfuer[0,:];  sy = esfuer[1,:];  txy = esfuer[2,:]
ex = deform[0,:];  ey = deform[1,:];  gxy = deform[2,:]
ez = -(nue/Ee)*(sx+sy)          # deformaciones ez en tensión plana

# %% imprimo y grafico las deformaciones
tabla_exeyezgxy = pd.DataFrame(
   data=np.c_[ex, ey, ez, gxy],
   index=np.arange(nef)+1,
   columns=['ex', 'ey', 'ez', 'gxy [rad]'])
tabla_exeyezgxy.index.name = '# EF'

compartir_variables(xnod, LaG, cg, interpolar=True)
plot_esf_def(ex,  r'$\epsilon_x$')
plot_esf_def(ey,  r'$\epsilon_y$')
plot_esf_def(ez,  r'$\epsilon_z$')
plot_esf_def(gxy, r'$\gamma_{xy}$ [rad]')

# %% imprimo y grafico los esfuerzos
tabla_sxsytxy = pd.DataFrame(
   data=np.c_[sx, sy, txy],
   index=np.arange(nef)+1,
   columns=['sx [Pa]', 'sy [Pa]', 'txy [Pa]'])
tabla_sxsytxy.index.name = '# EF'

plot_esf_def(sx,  r'$\sigma_x$ [Pa]')
plot_esf_def(sy,  r'$\sigma_y$ [Pa]')
plot_esf_def(txy, r'$\tau_{xy}$ [Pa]')

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

plot_esf_def(s1,   r'$\sigma_1$ [Pa]',   [ ang ])
plot_esf_def(s2,   r'$\sigma_2$ [Pa]',   [ ang+np.pi/2 ] )
plot_esf_def(tmax, r'$\tau_{max}$ [Pa]', [ ang-np.pi/4, ang+np.pi/4 ])
plot_esf_def(sv,   r'$\sigma_{VM}$ [Pa]')

# %%Se escriben los resultados a una hoja de MS EXCEL
archivo_resultados = f"resultados_{nombre_archivo}.xlsx"
writer = pd.ExcelWriter(archivo_resultados, engine='xlsxwriter')

# se escribe cada DataFrame a una hoja diferente
tabla_afq.to_excel(       writer, sheet_name='afq')
tabla_exeyezgxy.to_excel( writer, sheet_name='exeyezgxy')
tabla_sxsytxy.to_excel(   writer, sheet_name='sxsytxy')
tabla_s1s2tmaxsv.to_excel(writer, sheet_name='s1s2tmaxsv')

# Se cierra y graba el archivo de MS EXCEL
writer.save()
print(f'Cálculo finalizado. En "{archivo_resultados}" se guardaron los resultados.')

# %% Se genera un archivo .VTK para visualizar en Paraview
# Instale meshio (https://github.com/nschloe/meshio) con:
# pip install meshio[all] --user

import meshio
meshio.write_points_cells(
    f"resultados_{nombre_archivo}.vtk",
    points=xnod,
    cells={"triangle": LaG },
    point_data = {
        'uv'         : a.reshape((nno,2)),
        'reacciones' : q.reshape((nno,2))
        },
    cell_data = {
        "triangle" : 
        {
            'ex':ex, 'ey':ey, 'ez':ez,     'gxy':gxy,
            'sx':sx, 'sy':sy, 'txy':txy,
            's1':s1, 's2':s2, 'tmax':tmax, 'sv':sv,
            'n1':np.c_[np.cos(ang),           np.sin(ang)          ],
            'n2':np.c_[np.cos(ang + np.pi/2), np.sin(ang + np.pi/2)],
            "material" : mat
        }
    }
    # field_data=field_data
)

# %% bye, bye!
