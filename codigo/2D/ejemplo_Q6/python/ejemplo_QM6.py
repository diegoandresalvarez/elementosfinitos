# -*- coding: utf-8 -*-

#%%
'''
-------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para TENSION PLANA usando elementos
      rectangulares de 4 nodos con modos incompatibles
-------------------------------------------------------------------------------

DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramientos, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from funciones import t2ft_R4, compartir_variables, plot_esf_def

# %% constantes que ayudarán en la lectura del código
X, Y = 0, 1
NL1, NL2, NL3, NL4, NL5, NL6 = range(6)
r_ = [ 0, 1, 2, 3, 4, 5, 6, 7 ]   # GDL a retener en condensación nodal
e_ = [ 8, 9, 10, 11 ]             # GDL a eliminar en condensación nodal

g = 9.81 # [m/s²]   aceleración de la gravedad

# %% seleccione la malla a emplear:
#nombre_archivo = 'malla1'    # EJEMPLO CLASE
nombre_archivo = 'malla1_no_estructurada'    # EJEMPLO CLASE
df = pd.read_excel(f"{nombre_archivo}.xlsx", sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=0 o Y=1
xnod = df['xnod'][['x', 'y']].to_numpy()
nno  = xnod.shape[0]    # número de nodos (número de filas de la matriz xnod)

# %% definición de los grados de libertad
ngdl = 2*nno            # número de grados de libertad por nodo = [X, Y]
gdl  = np.reshape(np.arange(ngdl), (nno, 2)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG = df['LaG_mat'][['NL1', 'NL2', 'NL3', 'NL4']].to_numpy() - 1
nef = LaG.shape[0]      # número de EFs (número de filas de la matriz LaG)

# %% definición de los materiales
mat  = df['LaG_mat']['material'].to_numpy() - 1
Ee   = df['prop_mat']['E'].to_numpy()       # [Pa]     módulo de elasticidad del sólido
nue  = df['prop_mat']['nu'].to_numpy()      # [-]      coeficiente de Poisson
rhoe = df['prop_mat']['rho'].to_numpy()     # [kg/m³]  densidad
te   = df['prop_mat']['espesor'].to_numpy() # [m]      espesor
nmat = Ee.shape[0]                          # número de materiales

# %% relación de cargas puntuales
cp  = df['carga_punt']
ncp = cp.shape[0]       # número de cargas puntuales
f   = np.zeros(ngdl)    # vector de fuerzas nodales equivalentes global
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

# %% Se dibuja la malla de elementos finitos
cg = np.zeros((nef,2))  # almacena el centro de gravedad de los EF
plt.figure()
for e in range(nef):
    # se dibujan las aristas
    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
    plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], 'b')
    # se calcula la posición del centro de gravedad
    cg[e] = np.mean(xnod[LaG[e]], axis = 0)
    # y se reporta el número del elemento actual
    plt.text(cg[e,X], cg[e,Y], f'{e+1}', horizontalalignment='center',
                                         verticalalignment='center',  color='b')

# en todos los nodos se dibuja un marcador y se reporta su numeración
plt.plot(xnod[:, X], xnod[:, Y], 'r*')
for i in range(nno):
    plt.text(xnod[i, X], xnod[i, Y], f'{i+1}', color = 'r')

plt.gca().set_aspect('equal', adjustable = 'box')
plt.tight_layout()
plt.title('Malla de elementos finitos')
plt.show()

#%% Funciones de forma y sus derivadas del elemento rectangular de 4 nodos (N1 a
#   N4) y las funciones de forma de los modos incompatibles (N5 y N6)
Nforma = lambda xi,eta: np.array(
                        [ ((eta - 1)*(xi - 1))/4,    # N1
                         -((eta - 1)*(xi + 1))/4,    # N2
                          ((eta + 1)*(xi + 1))/4,    # N3
                         -((eta + 1)*(xi - 1))/4, ]) # N4
#                                      1 - xi**2,    # N5
#                                     1 - eta**2  ]) # N6

# derivadas de las funciones de forma con respecto a xi
dN_dxi = lambda xi,eta: np.array(
                        [  eta/4 - 1/4,         # dN1_dxi
                           1/4 - eta/4,         # dN2_dxi
                           eta/4 + 1/4,         # dN3_dxi
                         - eta/4 - 1/4,         # dN4_dxi
                                - 2*xi,         # dN5_dxi
                                     0     ])   # dN6_dxi

# derivadas de N con respecto a eta
dN_deta = lambda xi,eta: np.array(
                         [  xi/4 - 1/4,         # dN1_deta
                           -xi/4 - 1/4,         # dN2_deta
                            xi/4 + 1/4,         # dN3_deta
                            1/4 - xi/4,         # dN4_deta
                                     0,         # dN5_deta
                               - 2*eta      ])  # dN6_deta

#%% Cuadratura de Gauss-Legendre
# NOTA: se asumirá aquí el mismo orden de la cuadratura tanto en la dirección
#       de xi como en la dirección de eta
n_gl       = 2                       # orden de la cuadratura de Gauss-Legendre
x_gl, w_gl = np.polynomial.legendre.leggauss(n_gl)

# %% Ensamblaje la matriz de rigidez global y el vector de fuerzas másicas
#    nodales equivalentes global

# se inicializan la matriz de rigidez global y los espacios en memoria que
#  almacenarán las matrices de forma y de deformación
K       = np.zeros((ngdl, ngdl))          # matriz de rigidez global
inv_Kee = np.zeros((nef, 4, 4))
Ker     = np.zeros((nef, 4, 8))

N = np.empty((nef,n_gl,n_gl,2,2*4)) # matriz de forma en cada punto de GL
B = np.empty((nef,n_gl,n_gl,3,2*6)) # matriz de deformaciones en cada punto de GL
idx = nef * [None]                  # indices asociados a los gdl del EF e

# matriz constitutiva del elemento para TENSION PLANA
De = nmat * [ None ]
be = nmat * [ None ]
for i in range(nmat):
    De[i] = np.array([[Ee[i]/(1 - nue[i]**2),        Ee[i]*nue[i]/(1 - nue[i]**2), 0                     ],
                      [Ee[i]*nue[i]/(1 - nue[i]**2), Ee[i]/(1 - nue[i]**2),        0                     ],
                      [0,                            0,                            Ee[i]/(2*(1 + nue[i]))]])
    be[i] = np.array([0, -rhoe[i]*g])  # [kgf/m³] vector de fuerzas másicas

# para cada elemento finito en la malla:
for e in range(nef):
    # se calculan con el siguiente ciclo las matrices de rigidez y el vector de
    # fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
    Ke16 = np.zeros((12, 12))
    fe   = np.zeros(8)
    det_Je = np.empty((n_gl, n_gl))     # matriz para almacenar los jacobianos    
    
    # se llaman las coordenadas nodales del elemento para calcular las
    # derivadas de la función de transformación
    xe, ye = xnod[LaG[e], X], xnod[LaG[e], Y]
    
    # se calcula el área del cuadrilatero (del EF)
    # Ver https://en.wikipedia.org/wiki/Shoelace_formula
    # Ae = 0.5*(x1y2 + x2y3 + x3y4 + x4y1 - x2y1 - x3y2 - x4y3 - x1y4)
    Ae = 0.5*(xe[0]*ye[1] + xe[1]*ye[2] + xe[2]*ye[3] + xe[3]*ye[0] 
            - xe[1]*ye[0] - xe[2]*ye[1] - xe[3]*ye[2] - xe[0]*ye[3])


    for p in range(n_gl):
        for q in range(n_gl):
            # en cada punto de la cuadratura de Gauss-Legendre se evalúan las
            # funciones de forma y sus derivadas
            xi_gl, eta_gl = x_gl[p], x_gl[q]

            NNforma  = Nforma (xi_gl, eta_gl)
            ddN_dxi  = dN_dxi (xi_gl, eta_gl)
            ddN_deta = dN_deta(xi_gl, eta_gl)

            # se calcula el Jacobiano con N1 a N4
            dx_dxi  = np.sum(ddN_dxi[:4]  * xe);    dy_dxi  = np.sum(ddN_dxi[:4]  * ye)
            dx_deta = np.sum(ddN_deta[:4] * xe);    dy_deta = np.sum(ddN_deta[:4] * ye)

            # con ellas se ensambla la matriz Jacobiana del elemento y se
            # calcula su determinante
            Je = np.array([[dx_dxi,  dy_dxi ],
                           [dx_deta, dy_deta]])
            det_Je[p, q] = np.linalg.det(Je)

            # las matrices de forma y de deformación se evalúan y se ensamblan
            # en el punto de Gauss
            Npq = np.empty((2, 2*4))
            Bpq = np.empty((3, 2*6))
            for i in range(4):
                Npq[:,[2*i, 2*i+1]] = np.array([[NNforma[i], 0         ],
                                                [0,          NNforma[i]]])

            for i in range(6):
                dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je[p,q]
                dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je[p,q]
                Bpq[:,[2*i, 2*i+1]] = np.array([[dNi_dx, 0     ],
                                                [0,      dNi_dy],
                                                [dNi_dy, dNi_dx]])
            N[e,p,q] = Npq
            B[e,p,q] = Bpq

    # se calcula la matriz de correccion B56c    
    B56c = np.zeros((3,4))
    for p in range(n_gl):
        for q in range(n_gl):            
            B56c += -(1/Ae) * B[e,p,q][:,8:] * det_Je[p,q]*w_gl[p]*w_gl[q]

    for p in range(n_gl):
        for q in range(n_gl):
            # y se agrega la matriz de correccion B56c a B56
            B[e,p,q][:,8:] += B56c
            
            # se ensamblan la matriz de rigidez del elemento y el vector de
            # fuerzas nodales equivalentes del elemento
            Ke16 += B[e,p,q].T @ De[mat[e]] @ B[e,p,q] * det_Je[p,q]*te[mat[e]]*w_gl[p]*w_gl[q]
            fe   += N[e,p,q].T @ be[mat[e]]            * det_Je[p,q]*te[mat[e]]*w_gl[p]*w_gl[q]

    # se determina si hay puntos con jacobiano negativo, en caso tal se termina
    # el programa y se reporta
    if np.any(det_Je <= 0):
        raise Exception(f'Hay puntos con det_Je negativo en el elemento {e+1}')

    # se condensan los GDL jerárquicos u5, v5, u6, v6
    Krr = Ke16[np.ix_(r_,r_)];  Ker[e]     = Ke16[np.ix_(e_,r_)]
    Kre = Ke16[np.ix_(r_,e_)];  inv_Kee[e] = np.linalg.inv(Ke16[np.ix_(e_,e_)])
    Ke = Krr - Kre @ inv_Kee[e] @ Ker[e]

    # y se añaden la matriz de rigidez del elemento y el vector de fuerzas
    # nodales del elemento a sus respectivos arreglos de la estructura
    idx[e] = gdl[LaG[e]].flatten() # se obtienen los grados de libertad
    K[np.ix_(idx[e], idx[e])] += Ke
    f[np.ix_(idx[e])]         += fe

# %% Muestro la configuración de la matriz K (K es rala)
plt.figure()
plt.spy(K)
plt.title('Los puntos representan los elementos diferentes de cero')
plt.show()

#%% Cálculo de las cargas nodales equivalentes de las cargas distribuidas:
cd   = df['carga_distr']
nlcd = cd.shape[0]     # número de lados con carga distribuida
ft   = np.zeros(ngdl)  # fuerzas nodales equivalentes de cargas superficiales

# por cada lado cargado se obtienen las fuerzas nodales equivalentes en los
# nodos y se añaden al vector de fuerzas superficiales
for i in range(nlcd):
   e     = cd['elemento'][i] - 1
   lado  = cd['lado'][i]
   carga = cd[['tix', 'tiy', 'tjx', 'tjy']].loc[i].to_numpy()
   fte = t2ft_R4(xnod[LaG[e,:],:], lado, carga, te[mat[e]])
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
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # separo la memoria
a[c] = ac;          a[d] = ad          # desplazamientos
q[c] = qd         # q[d] = qc = 0      # fuerzas nodales de equilibrio

# %% Dibujo la malla de elementos finitos y las deformada de esta
delta  = np.reshape(a, (nno,2))
escala = 50000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

plt.figure()
for e in range(nef):
   nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
   plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], 'r',
                        label='Posición original'  if e == 0 else "", lw=0.5)
   plt.plot(xdef[nod_ef, X], xdef[nod_ef, Y], 'b',
                        label='Posición deformada' if e == 0 else "")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.title(f'Deformada escalada {escala} veces')
plt.tight_layout()
plt.show()

#%% Deformaciones y los esfuerzos en los puntos de Gauss
deform = np.zeros((nef,n_gl,n_gl,3)) # deformaciones en cada punto de GL
esfuer = np.zeros((nef,n_gl,n_gl,3)) # esfuerzos en cada punto de GL

for e in range(nef):
    ar = a[idx[e]]    # desplazamientos nodales del elemento e
    ae = np.r_[ ar,
               -inv_Kee[e] @ Ker[e] @ ar]
    for pp in range(n_gl):
        for qq in range(n_gl):
            deform[e,pp,qq] = B[e,pp,qq] @ ae              # calculo las deformaciones
            esfuer[e,pp,qq] = De[mat[e]] @ deform[e,pp,qq] # calculo los esfuerzos
    '''
    ae = a[idx[e]]    # desplazamientos nodales del elemento e
    for pp in range(n_gl):
        for qq in range(n_gl):
            deform[e,pp,qq] = B[e,pp,qq][:,:8] @ ae        # calculo las deformaciones
            esfuer[e,pp,qq] = De[mat[e]] @ deform[e,pp,qq] # calculo los esfuerzos
    '''

#%% Esfuerzos y deformaciones en los nodos:
num_elem_ady = np.zeros(nno)
sx  = np.zeros(nno);    txy = np.zeros(nno);    ex  = np.zeros(nno)
sy  = np.zeros(nno);    txz = np.zeros(nno);    ey  = np.zeros(nno)
sz  = np.zeros(nno);    tyz = np.zeros(nno);    gxy = np.zeros(nno)

# matriz de extrapolación
A = np.array([
    [  3**(1/2)/2 + 1,             -1/2,             -1/2,   1 - 3**(1/2)/2],
    [            -1/2,   1 - 3**(1/2)/2,   3**(1/2)/2 + 1,             -1/2],
    [  1 - 3**(1/2)/2,             -1/2,             -1/2,   3**(1/2)/2 + 1],
    [            -1/2,   3**(1/2)/2 + 1,   1 - 3**(1/2)/2,             -1/2]])

# se hace la extrapolación de los esfuerzos y las deformaciones en cada elemento
# a partir de las lecturas en los puntos de Gauss
for e in range(nef):
    #sx[LaG[e]]  += A @ np.array([esfuer[e,0,0,0],   # I   = (p=0, q=0)
    #                             esfuer[e,0,1,0],   # II  = (p=0, q=1)
    #                             esfuer[e,1,0,0],   # III = (p=1, q=0)
    #                             esfuer[e,1,1,0]])  # IV  = (p=1, q=1)
    sx [LaG[e]] += A @ esfuer[e,:,:,0].ravel()
    sy [LaG[e]] += A @ esfuer[e,:,:,1].ravel()
    txy[LaG[e]] += A @ esfuer[e,:,:,2].ravel()
    ex [LaG[e]] += A @ deform[e,:,:,0].ravel()
    ey [LaG[e]] += A @ deform[e,:,:,1].ravel()
    gxy[LaG[e]] += A @ deform[e,:,:,2].ravel()

    # se lleva un conteo de los elementos adyacentes a un nodo
    num_elem_ady[LaG[e]] += 1

# en todos los nodos se promedia los esfuerzos y las deformaciones de los
# elementos, se alisa la malla de resultados
sx  /= num_elem_ady;   ex  /= num_elem_ady
sy  /= num_elem_ady;   ey  /= num_elem_ady
txy /= num_elem_ady;   gxy /= num_elem_ady

# se calculan las deformaciones ez en tension plana
ez = -(nue/Ee)*(sx + sy)

# %% Se calculan y grafican para cada elemento los esfuerzos principales y
#    sus direcciones
# NOTA: esto solo es valido para el caso de TENSION PLANA).
# En caso de DEFORMACIÓN PLANA se deben calcular los valores y vectores
# propios de la matriz de tensiones de Cauchy
#   [dirppales{e}, esfppales{e}] = eigh([sx  txy 0    % matriz de esfuerzos
#                                        txy sy  0    % de Cauchy
#                                        0   0   0]);
s1   = (sx+sy)/2 + np.sqrt(((sx-sy)/2)**2 + txy**2) # esfuerzo normal máximo
s2   = (sx+sy)/2 - np.sqrt(((sx-sy)/2)**2 + txy**2) # esfuerzo normal mínimo
tmax = (s1 - s2)/2                                  # esfuerzo cortante máximo
ang  = 0.5*np.arctan2(2*txy, sx-sy)                 # ángulo de inclinación de s1

# %% Calculo de los esfuerzos de von Mises
s3   = np.zeros(nno)
sv   = np.sqrt(((s1-s2)**2 + (s2-s3)**2 + (s1-s3)**2)/2)

# %% Gráfica del post-proceso:
# las matrices xnod y LaG se vuelven variables globales por facilidad
compartir_variables(xnod, LaG)

# deformaciones
plot_esf_def(ex,   r'$\epsilon_x$')
plot_esf_def(ey,   r'$\epsilon_y$')
plot_esf_def(ez,   r'$\epsilon_z$')
plot_esf_def(gxy,  r'$\gamma_{xy}$ [rad]')

# esfuerzos
plot_esf_def(sx,   r'$\sigma_x$ [Pa]')
plot_esf_def(sy,   r'$\sigma_y$ [Pa]')
plot_esf_def(txy,  r'$\tau_{xy}$ [Pa]')

# esfuerzos principales con sus orientaciones
plot_esf_def(s1,   r'$\sigma_1$ [Pa]',     ang                       )
plot_esf_def(s2,   r'$\sigma_2$ [Pa]',     ang+np.pi/2               )
plot_esf_def(tmax, r'$\tau_{máx}$ [Pa]', [ ang-np.pi/4, ang+np.pi/4 ])

# esfuerzos de von Mises
plot_esf_def(sv,   r'$\sigma_{VM}$ [Pa]')

# %% Reporte de los resultados:

# se crean tablas para reportar los resultados nodales de: desplazamientos (a),
# fuerzas nodales equivalentes (f) y fuerzas nodales de equilibrio (q)
tabla_afq = pd.DataFrame(
    data=np.c_[a.reshape((nno,2)), f.reshape((nno,2)), q.reshape((nno,2))],
    index=np.arange(nno)+1,
    columns=['ux [m]', 'uy [m]', 'fx [N]', 'fy [N]', 'qx [N]', 'qy [N]'])
tabla_afq.index.name = '# nodo'

# deformaciones
tabla_def = pd.DataFrame(data    = np.c_[ex, ey, ez, gxy],
                         index   = np.arange(nno) + 1,
                         columns = ['ex', 'ey', 'ez', 'gxy [rad]'])
tabla_def.index.name = '# nodo'

# esfuerzos
tabla_esf = pd.DataFrame(data    = np.c_[sx, sy, txy],
                         index   = np.arange(nno) + 1,
                         columns = ['sx [Pa]', 'sy [Pa]', 'txy [Pa]'])
tabla_esf.index.name = '# nodo'

# esfuerzos principales y de von Misses:
tabla_epv = pd.DataFrame(
       data    = np.c_[s1, s2, tmax, sv, ang],
       index   = np.arange(nno) + 1,
       columns = ['s1 [Pa]', 's2 [Pa]', 'tmax [Pa]', 'sv [Pa]', 'theta [rad]'])
tabla_epv.index.name = '# nodo'

# se crea un archivo de MS EXCEL
archivo_resultados = f"resultados_{nombre_archivo}.xlsx"
writer = pd.ExcelWriter(archivo_resultados, engine = 'xlsxwriter')

# cada tabla hecha previamente es guardada en una hoja del archivo de Excel
tabla_afq.to_excel(writer, sheet_name = 'afq')
tabla_def.to_excel(writer, sheet_name = 'exeyezgxy')
tabla_esf.to_excel(writer, sheet_name = 'sxsytxy')
tabla_epv.to_excel(writer, sheet_name = 's1s2tmaxsv')
writer.save()

print(f'Cálculo finalizado. En "{archivo_resultados}" se guardaron los resultados.')

# %% Se genera un archivo .VTK para visualizar en Paraview
# Instale meshio (https://github.com/nschloe/meshio) con:
# pip install meshio[all] --user

import meshio
meshio.write_points_cells(
    f"resultados_{nombre_archivo}.vtk",
    points=xnod,
    cells={"quad": LaG[:,[0,1,2,3]] },
    point_data = {
        'ex':ex, 'ey':ey, 'ez'  :ez,   'gxy':gxy,
        'sx':sx, 'sy':sy, 'txy' :txy,
        's1':s1, 's2':s2, 'tmax':tmax, 'sv':sv,
        'uv':a.reshape((nno,2)),
        'reacciones' : q.reshape((nno,2)),
        'n1':np.c_[np.cos(ang),           np.sin(ang)          ],
        'n2':np.c_[np.cos(ang + np.pi/2), np.sin(ang + np.pi/2)]
        },
    cell_data = {
        "quad" : {"material" : mat}
    },
)

# %%bye, bye!
