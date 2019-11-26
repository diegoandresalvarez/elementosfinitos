# -*- coding: utf-8 -*-

#%%
'''
--------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para EFs hexaédricos de 20 nodos (H20)
--------------------------------------------------------------------------------

DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramientos, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse 
from scipy.sparse.linalg import spsolve
from funciones import gausslegendre_quad_hexa, matriz_extrapolacion_esfuerzos_H20
from funciones import t2ft_H20, N_H20, dN_dxi_H20, dN_deta_H20, dN_dzeta_H20

# %% constantes que ayudarán en la lectura del código
X, Y, Z = 0, 1, 2
g = 9.81 # [m/s²]   aceleración de la gravedad

# %% seleccione la malla a emplear:
nombre_archivo = 'malla_H20_viga'
# nombre_archivo = 'malla_H20_conexion'
df = pd.read_excel(f"{nombre_archivo}.xlsx", sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=0, Y=1, Z=2
xnod = df['xnod'][['x', 'y', 'z']].to_numpy()
nno  = xnod.shape[0]    # número de nodos (número de filas de la matriz xnod)

# %% definición de los grados de libertad
ngdl = 3*nno            # número de grados de libertad por nodo = [X, Y, Z]
gdl  = np.reshape(np.arange(ngdl), (nno, 3)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG  = df['LaG_mat'].loc[:, 'NL1':'NL20'].to_numpy() - 1
# se reorganizan las columnas con mi numeración local de los nodos. Esto es para 
# poder importar la malla, la cual fue generada en GiD
idx_LaG = np.array([1, 9, 2, 10, 3, 11, 4, 12, 13, 14, 15, 16, 5, 17, 6, 18, 7, 19, 8, 20]) - 1
LaG = LaG[:, idx_LaG]

# se carga el número del material
mat  = df['LaG_mat']['material'].to_numpy() - 1
nef  = LaG.shape[0]    # número de EFs (número de filas de la matriz LaG)
nnpe = LaG.shape[1]    # número de nodos por EF (=20)

if nnpe != 20:
   raise Exception('Este código SOLO es sirve para EFs hexaédricos de 20 nodos (H20)')

# %% material
Ee   = df['prop_mat']['E'].to_numpy()   # [Pa]     módulo de elasticidad
nue  = df['prop_mat']['nu'].to_numpy()  # [-]      coeficiente de Poisson
rhoe = df['prop_mat']['rho'].to_numpy() # [kg/m³]  densidad
nmat = Ee.shape[0]                      # número de materiales

# %% relación de cargas puntuales
cp  = df['carga_punt']
ncp = cp.shape[0]       # número de cargas puntuales
f   = np.zeros(ngdl)    # vector de fuerzas nodales equivalentes global
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

#%% Funciones de forma (serendípitas) y sus derivadas del elemento rectangular
#   de 8 nodos:
Nforma   = N_H20
dN_dxi   = dN_dxi_H20
dN_deta  = dN_deta_H20
dN_dzeta = dN_dzeta_H20

#%% Cuadratura de Gauss-Legendre
# NOTA: se asumirá aquí el mismo orden de la cuadratura tanto en la dirección
#       de xi como en la dirección de eta
x_gl, w_gl = gausslegendre_quad_hexa(2) # 2x2x2
n_gl = len(w_gl)

# %% Ensamblaje la matriz de rigidez global y el vector de fuerzas másicas
#    nodales equivalentes global

# se inicializan la matriz de rigidez global y los espacios en memoria que
#  almacenarán las matrices de forma y de deformación
#K= np.zeros((ngdl, ngdl))        # matriz de rigidez global
K = sparse.coo_matrix((ngdl, ngdl)) # matriz de rigidez global como RALA (sparse)

N = np.empty((nef,n_gl,3,3*nnpe)) # matriz de forma en cada punto de GL
B = np.empty((nef,n_gl,6,3*nnpe)) # matriz de deformaciones en cada punto de GL
idx = nef * [None]                # indices asociados a los gdl del EF e

# matriz constitutiva del elemento para TENSION PLANA
De = nmat * [ None ]
be = nmat * [ None ]
for i in range(nmat):
    d1 = (1-nue[i])/(1-2*nue[i])
    d2 = nue[i]/(1-2*nue[i])
    De[i] = Ee[i]/(1+nue[i])*np.array([[d1, d2, d2, 0,   0,   0  ],
                                       [d2, d1, d2, 0,   0,   0  ],
                                       [d2, d2, d1, 0,   0,   0  ],
                                       [0,  0,  0,  1/2, 0,   0  ],
                                       [0,  0,  0,  0,   1/2, 0  ],
                                       [0,  0,  0,  0,   0,   1/2]])
    be[i] = np.array([0, -rhoe[i]*g, 0])  # [kgf/m³] vector de fuerzas másicas

# para cada elemento finito en la malla:
for e in range(nef):
    # se calculan con el siguiente ciclo las matrices de rigidez y el vector de
    # fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
    Ke = np.zeros((3*nnpe, 3*nnpe))
    fe = np.zeros(3*nnpe)
    det_Je = np.empty(n_gl)      # matriz para almacenar los jacobianos

    for p in range(n_gl):
        # en cada punto de la cuadratura de Gauss-Legendre se evalúan las
        # funciones de forma y sus derivadas
        xi_gl, eta_gl, zeta_gl = x_gl[p,:]

        NNforma   = Nforma  (xi_gl, eta_gl, zeta_gl)
        ddN_dxi   = dN_dxi  (xi_gl, eta_gl, zeta_gl)
        ddN_deta  = dN_deta (xi_gl, eta_gl, zeta_gl)
        ddN_dzeta = dN_dzeta(xi_gl, eta_gl, zeta_gl)

        # se llaman las coordenadas nodales del elemento para calcular las
        # derivadas de la función de transformación
        xe, ye, ze = xnod[LaG[e], X], xnod[LaG[e], Y], xnod[LaG[e], Z]

        dx_dxi   = np.sum(ddN_dxi  * xe);    dy_dxi   = np.sum(ddN_dxi  * ye)
        dx_deta  = np.sum(ddN_deta * xe);    dy_deta  = np.sum(ddN_deta * ye)        
        dx_dzeta = np.sum(ddN_dzeta* xe);    dy_dzeta = np.sum(ddN_dzeta* ye)                

        dz_dxi   = np.sum(ddN_dxi  * ze)
        dz_deta  = np.sum(ddN_deta * ze)
        dz_dzeta = np.sum(ddN_dzeta* ze)

        # con ellas se ensambla la matriz Jacobiana del elemento y se
        # calcula su determinante
        Je = np.array([[dx_dxi,   dy_dxi,   dz_dxi  ],
                       [dx_deta,  dy_deta,  dz_deta ],
                       [dx_dzeta, dy_dzeta, dz_dzeta]])
        det_Je[p] = np.linalg.det(Je)

        # las matrices de forma y de deformación se evalúan y se ensamblan
        # en el punto de Gauss
        Np = np.empty((3, 3*nnpe))
        Bp = np.empty((6, 3*nnpe))
        for i in range(nnpe):
            Np[:,[3*i, 3*i+1, 3*i+2]] = np.array(
                [[NNforma[i], 0,          0         ],
                 [0,          NNforma[i], 0         ],
                 [0,          0,          NNforma[i]]])

            dNi_dx, dNi_dy, dNi_dz = np.linalg.solve(Je, 
                            np.array([ddN_dxi[i],  ddN_deta[i],  ddN_dzeta[i]]))

            Bp[:,[3*i, 3*i+1, 3*i+2]] = np.array(
                [[dNi_dx, 0,      0     ],
                 [0,      dNi_dy, 0     ],
                 [0,      0,      dNi_dz],       
                 [dNi_dy, dNi_dx, 0     ],
                 [dNi_dz, 0,      dNi_dx],
                 [0,      dNi_dz, dNi_dy ]])
        N[e,p] = Np
        B[e,p] = Bp

        # se ensamblan la matriz de rigidez del elemento y el vector de
        # fuerzas nodales equivalentes del elemento
        Ke += Bp.T @ De[mat[e]] @ Bp * det_Je[p]*w_gl[p]
        fe += Np.T @ be[mat[e]]      * det_Je[p]*w_gl[p]

    # se determina si hay puntos con jacobiano negativo, en caso tal se termina
    # el programa y se reporta
    if np.any(det_Je <= 0):
        raise Exception(f'Hay puntos con det_Je negativo en el elemento {e+1}')

    # y se añaden la matriz de rigidez del elemento y el vector de fuerzas
    # nodales del elemento a sus respectivos arreglos de la estructura
    idx[e] = gdl[LaG[e]].flatten() # se obtienen los grados de libertad
    # K[np.ix_(idx[e], idx[e])] += Ke
    IDX = np.array([(i,j) for i in idx[e] for j in idx[e]])
    K += sparse.coo_matrix((Ke.flat, (IDX[:,0],IDX[:,1])), shape=(ngdl,ngdl))

    f[np.ix_(idx[e])]         += fe

'''
# %% Muestro la configuración de la matriz K (K es rala)
plt.figure()
plt.spy(K)
plt.title('Los puntos representan los elementos diferentes de cero')
plt.show()
'''

#%% Cálculo de las cargas nodales equivalentes de las cargas distribuidas:
cd   = df['carga_distr']
nlcd = cd.shape[0]     # número de lados con carga distribuida
ft   = np.zeros(ngdl)  # fuerzas nodales equivalentes de cargas superficiales

'''
# por cada lado cargado se obtienen las fuerzas nodales equivalentes en los
# nodos y se añaden al vector de fuerzas superficiales
for i in range(nlcd):
   e     = cd['elemento'][i] - 1
   lado  = cd['lado'][i]
   carga = cd[['tix', 'tiy', 'tjx', 'tjy', 'tkx', 'tky']].loc[i].to_numpy()
   fte   = t2ft_H20(xnod[LaG[e,:],:], lado, carga, te)

   ft[np.ix_(idx[e])] += fte
'''

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
#ad= np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
ad = spsolve(Kdd, fc - Kdc@ac)

qd = Kcc@ac + Kcd@ad - fd              # fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # separo la memoria
a[c] = ac;          a[d] = ad          # desplazamientos
q[c] = qd         # q[d] = qc = 0      # fuerzas nodales de equilibrio

#%% Deformaciones y los esfuerzos en los puntos de Gauss
deform = np.zeros((nef,n_gl,6)) # deformaciones en cada punto de GL
esfuer = np.zeros((nef,n_gl,6)) # esfuerzos en cada punto de GL

for e in range(nef):
    ae = a[idx[e]]    # desplazamientos nodales del elemento e
    for p in range(n_gl):
        deform[e,p] = B[e,p] @ ae              # calculo las deformaciones
        esfuer[e,p] = De[mat[e]] @ deform[e,p] # calculo los esfuerzos

#%% Esfuerzos y deformaciones en los nodos:
num_elem_ady = np.zeros(nno)
sx  = np.zeros(nno);        ex  = np.zeros(nno)
sy  = np.zeros(nno);        ey  = np.zeros(nno)
sz  = np.zeros(nno);        ez  = np.zeros(nno)
txy = np.zeros(nno);        gxy = np.zeros(nno)
txz = np.zeros(nno);        gxz = np.zeros(nno)
tyz = np.zeros(nno);        gyz = np.zeros(nno)

# matriz de extrapolación
A = matriz_extrapolacion_esfuerzos_H20()

# se hace la extrapolación de los esfuerzos y las deformaciones en cada elemento
# a partir de las lecturas en los puntos de Gauss
for e in range(nef):
    sx [LaG[e]] += A @ esfuer[e,:,0].ravel()
    sy [LaG[e]] += A @ esfuer[e,:,1].ravel()
    sz [LaG[e]] += A @ esfuer[e,:,2].ravel()
    txy[LaG[e]] += A @ esfuer[e,:,3].ravel()
    txz[LaG[e]] += A @ esfuer[e,:,4].ravel()    
    tyz[LaG[e]] += A @ esfuer[e,:,5].ravel()        
    ex [LaG[e]] += A @ deform[e,:,0].ravel()
    ey [LaG[e]] += A @ deform[e,:,1].ravel()
    ez [LaG[e]] += A @ deform[e,:,2].ravel()    
    gxy[LaG[e]] += A @ deform[e,:,3].ravel()
    gxz[LaG[e]] += A @ deform[e,:,4].ravel()
    gyz[LaG[e]] += A @ deform[e,:,5].ravel()    

    # se lleva un conteo de los elementos adyacentes a un nodo
    num_elem_ady[LaG[e]] += 1

# en todos los nodos se promedia los esfuerzos y las deformaciones de los
# elementos, se alisa la malla de resultados
sx  /= num_elem_ady;   ex  /= num_elem_ady
sy  /= num_elem_ady;   ey  /= num_elem_ady
sz  /= num_elem_ady;   ez  /= num_elem_ady
txy /= num_elem_ady;   gxy /= num_elem_ady
txz /= num_elem_ady;   gxz /= num_elem_ady
tyz /= num_elem_ady;   gyz /= num_elem_ady

# %% Se calculan para cada nodo los esfuerzos principales y sus direcciones
s1 = np.zeros(nno);  n1 = np.zeros((nno, 3))
s2 = np.zeros(nno);  n2 = np.zeros((nno, 3))
s3 = np.zeros(nno);  n3 = np.zeros((nno, 3))
for i in range(nno):
   esfppales, dirppales = np.linalg.eigh(
                             [[sx[i],   txy[i],  txz[i]],  # matriz de esfuerzos
                              [txy[i],  sy[i],   tyz[i]],  # de Cauchy
                              [txz[i],  tyz[i],  sz[i] ]])

   idx_esf = esfppales.argsort()[::-1] # ordene de mayor a menor
   s1[i], s2[i], s3[i] = esfppales[idx_esf]
   n1[i] = dirppales[:,idx_esf[0]]
   n2[i] = dirppales[:,idx_esf[1]]
   n3[i] = dirppales[:,idx_esf[2]]

# Esfuerzo cortante máximo
tmax = (s1-s3)/2                               # esfuerzo cortante máximo
   
# %% Calculo de los esfuerzos de von Mises
sv   = np.sqrt(((s1-s2)**2 + (s2-s3)**2 + (s1-s3)**2)/2)

# %% Reporte de los resultados:

# se crean tablas para reportar los resultados nodales de: desplazamientos (a),
# fuerzas nodales equivalentes (f) y fuerzas nodales de equilibrio (q)
tabla_afq = pd.DataFrame(
    data=np.c_[a.reshape((nno,3)), f.reshape((nno,3)), q.reshape((nno,3))],
    index=np.arange(nno)+1,
    columns=['u [m]',  'v [m]',  'w [m]',       
             'fx [N]', 'fy [N]', 'fz [N]',
             'qx [N]', 'qy [N]', 'qz [N]'])
tabla_afq.index.name = '# nodo'

# deformaciones
tabla_def = pd.DataFrame(data    = np.c_[ex, ey, ez, gxy, gxz, gyz],
                         index   = np.arange(nno) + 1,
                         columns = ['ex',        'ey',        'ez', 
                                    'gxy [rad]', 'gxz [rad]', 'gyz [rad]'])
tabla_def.index.name = '# nodo'

# esfuerzos
tabla_esf = pd.DataFrame(data    = np.c_[sx, sy, sz, txy, txz, tyz],
                         index   = np.arange(nno) + 1,
                         columns = ['sx [Pa]',  'sy [Pa]',  'sz [Pa]', 
                                    'txy [Pa]', 'txz [Pa]', 'tyz [Pa]'])
tabla_esf.index.name = '# nodo'

# esfuerzos principales y de von Misses:
tabla_epv = pd.DataFrame(
       data    = np.c_[s1, s2, s3, tmax, sv, n1, n2, n3],
       index   = np.arange(nno) + 1,
       columns = ['s1 [Pa]', 's2 [Pa]', 's3 [Pa]', 'tmax [Pa]', 'sv [Pa]',
                  'n1x', 'n1y', 'n1z',
                  'n2x', 'n2y', 'n2z', 
                  'n3x', 'n3y', 'n3z'])
tabla_epv.index.name = '# nodo'

# se crea un archivo de MS EXCEL
archivo_resultados = f"resultados/{nombre_archivo}.xlsx"
writer = pd.ExcelWriter(archivo_resultados, engine = 'xlsxwriter')

# cada tabla hecha previamente es guardada en una hoja del archivo de Excel
tabla_afq.to_excel(writer, sheet_name = 'afq')
tabla_def.to_excel(writer, sheet_name = 'deformaciones')
tabla_esf.to_excel(writer, sheet_name = 'esfuerzos')
tabla_epv.to_excel(writer, sheet_name = 'esf_ppales')
writer.save()

print(f'Cálculo finalizado. En "{archivo_resultados}" se guardaron los resultados.')

# %% Se genera un archivo .VTK para visualizar en Paraview
# Instale meshio (https://github.com/nschloe/meshio) con:
# pip install meshio[all] --user

# tenga en cuanta que VTK tiene una numeración diferente de los nodos que la 
# especificada en el libro de Oñate y que la usada en GiD
# https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

import meshio
meshio.write_points_cells(
    f"resultados/{nombre_archivo}.vtk",
    points=xnod,
    cells={"hexahedron20": LaG[:,np.array([3, 5, 7, 1, 15, 17, 19, 13, 4, 6, 8, 2, 16, 18, 20, 14, 10, 11, 12, 9])-1] },
    point_data = {
        'ex':ex, 'ey':ey, 'ez':ez, 'gxy':gxy,   'gxz':gxz, 'gyz':gyz,
        'sx':sx, 'sy':sy, 'sz':sz, 'txy':txy,   'txz':txz, 'tyz':tyz,
        's1':s1, 's2':s2, 's3':s3, 'tmax':tmax, 'sv':sv,
        'uvw' : a.reshape((nno,3)),
        'n1'  : n1, 
        'n2'  : n2, 
        'n3'  : n3
        }
    # cell_data=cell_data,
    # field_data=field_data
)

# %% Pasando los resultados a GiD
# Pasando los esfuerzos ya promediados:
# export_to_GiD('resultados/conexion_esf_nodos',xnod,LaG,a,q,[sx sy sz txy txz tyz])

# Pasando los puntos de Gauss [RECOMENDADO] !!!
# export_to_GiD('resultados/conexion_H20_esf_GP',xnod,LaG,a,q,esf);

# %%bye, bye!
