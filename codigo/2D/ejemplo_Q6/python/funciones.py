# -*- coding: utf-8 -*-

import numpy                as np
import matplotlib.pyplot    as plt
import warnings

# %% constantes que ayudaran en la lectura del codigo
X, Y = 0, 1
NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = range(8)

# %% variables globales que se heredarán del programa principal
xnod = None
LaG = None

# %%
def compartir_variables(xnod_, LaG_):
    '''
    Importa variables globales del programa principal a este módulo
    '''
    global xnod, LaG
    xnod, LaG = xnod_, LaG_
        
# %%
def t2ft_R4(xnod, lado, carga, espesor):
    '''Función que convierte las fuerzas superficiales aplicadas a un elemento
    finito rectangular de 4 nodos a sus correspondientes cargas nodales 
    equivalentes ft
    
    Recibe:
        xnod:  coordenadas nodales del elemento finito de 4 nodos
          xnod = [ x1e y1e
                   x2e y2e
                   x3e y3e
                   x4e y4e ]

        lado:  arista en la que se aplica la carga, puede tomar los siguientes
               valores: 12, 23, 34, 41

        carga: fuerza distribuida en los nodos
        
               [ t1x t1y t2x t2y ]; % si carga se aplica sobre lado 12
               [ t2x t2y t3x t3y ]; % si carga se aplica sobre lado 23
               [ t3x t3y t4x t4y ]; % si carga se aplica sobre lado 34
               [ t4x t4y t1x t1y ]; % si carga se aplica sobre lado 41
    
        espesor: espesor del elemento
    '''
    
    # se definen los indices de los lados
    if   lado == 12: idx = np.array([ 1, 2 ]) - 1
    elif lado == 23: idx = np.array([ 2, 3 ]) - 1
    elif lado == 34: idx = np.array([ 3, 4 ]) - 1
    elif lado == 41: idx = np.array([ 4, 1 ]) - 1
    else: 
        raise Exception('Únicamente se permiten los lados 12, 23, 34 o 41')

    nno = xnod.shape[0]
    if nno != 4:
        raise Exception('Solo para elementos rectangulares de 4 nodos')

    # parámetros para mejorar la lectura del código
    X, Y = 0, 1
   
    # se define el número de puntos de la cuadratura y se obtienen los puntos
    # de evaluación y los pesos de la cuadratura
    n_gl       = 5
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gl)
    
    # se definen las funciones de forma unidimensionales y sus derivadas
    NN      = lambda xi: np.array([ (1-xi)/2, (1+xi)/2 ])
    dNN_dxi = lambda xi: np.array([     -1/2,      1/2 ])
       
    # se calcula el vector de fuerzas distribuidas en los nodos
    te = np.zeros(2*nno)
    te[np.c_[2*idx, 2*idx + 1].ravel()] = carga
    
    # cálculo de la integral:
    suma   = np.zeros((2*nno, 2*nno))
    N      = np.zeros(nno)
    dN_dxi = np.zeros(nno)
    for p in range(n_gl):
        # se evalúan las funciones de forma
        N[idx] = NN(x_gl[p])
        matN = np.empty((2,2*nno))
        for i in range(nno):
            matN[:,[2*i, 2*i+1]] = np.array([[N[i], 0   ],
                                             [0,    N[i]]])

        # se calcula el jacobiano
        dN_dxi[idx] = dNN_dxi(x_gl[p])
        dx_dxi      = np.dot(dN_dxi, xnod[:,X])
        dy_dxi      = np.dot(dN_dxi, xnod[:,Y])
        ds_dxi      = np.hypot(dx_dxi, dy_dxi)

        # y se calcula la sumatoria
        suma += matN.T @ matN * ds_dxi*w_gl[p]

    # finalmente, se retorna el vector de fuerzas nodales equivalentes
    return espesor * (suma @ te)

#%%
def plot_esf_def(variable, titulo, angulo = None):
    '''
    Grafica variable para la malla de EFs especificada.

    Uso:
        variable: es la variable que se quiere graficar
        titulo:   título del gráfico
        angulo:   opcional para los esfuerzos principales s1 y s2 y tmax
    '''
    
    # Para propósitos de graficación el EF se divide en 2 triángulos así: 
    #     
    #                             7 -------6--------5
    #                             |       /|\       |
    #                             | EFT6 / | \ EFT3 |
    #                             |     /  |  \     |
    #                             |    /   |   \    |
    #                             |   /    |    \   |
    #                             |  /     |     \  |
    #                             | /      |      \ |
    #                             8/  EFT5 | EFT4  \4
    #                             |\       |       /|
    #                             | \      |      / |
    #                             |  \     |     /  |
    #                             |   \    |    /   |
    #                             |    \   |   /    |
    #                             |     \  |  /     |
    #                             | EFT1 \ | / EFT2 |
    #                             |       \|/       |
    #                             1--------2--------3
    
    

    #                             4--------3
    #                             |\       |
    #                             | \ EFT2 |
    #                             |  \     |
    #                             |   \    |
    #                             |    \   |
    #                             |     \  |
    #                             | EFT1 \ |
    #                             |       \|
    #                             1--------2
    
        
    nef = LaG.shape[0]    # número de elementos finitos en la malla original
    
    # se arma la matriz de correspondencia (LaG) de la nueva malla triangular
    LaG_t = np.zeros((2*nef, 3), dtype = int)
    for e in range(nef):
        LaG_t[2*e + 0, :] = LaG[e, [NL1, NL2, NL4]]
        LaG_t[2*e + 1, :] = LaG[e, [NL2, NL3, NL4]]
    
    # se inicializa el lienzo
    fig, ax = plt.subplots() 

    # se encuentra el máximo en valor absoluto para ajustar el colorbar()
    val_max = np.max(np.abs(variable))    

    # se grafica la malla de EFS, los colores en cada triángulo y las curvas 
    # de nivel
    for e in range(nef):
        # se dibujan las aristas
        nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
        plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.5, color = 'gray')

    im = ax.tripcolor(xnod[:, X], xnod[:, Y], LaG_t, variable, cmap = 'bwr',
                      shading = 'gouraud', vmin = -val_max, vmax = val_max)
    ax.tricontour(xnod[:, X], xnod[:, Y], LaG_t, variable, 20)
    
    # a veces sale un warning, simplemente porque no existe la curva 0
    warnings.filterwarnings("ignore")
    ax.tricontour(xnod[:, X], xnod[:, Y], LaG_t, variable, levels=[0], linewidths=3)
    warnings.filterwarnings("default")

    fig.colorbar(im, ax = ax, format = '%6.3g')
                
    # para los esfuerzos principales se grafican las líneas que indiquen las
    # direcciones de los esfuerzos en cada nodo de la malla
    if angulo is not None:
       if type(angulo) is np.ndarray: 
           angulo = [ angulo ]
       for ang in angulo:
           ax.quiver(xnod[:, X], xnod[:, Y], 
                variable*np.cos(ang), variable*np.sin(ang), 
                headwidth=0, headlength=0, headaxislength=0, pivot='middle')
    
    # se especifican los ejes y el título, y se colocan los ejes iguales
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')
    ax.set_title(titulo, fontsize=20)

    ax.set_aspect('equal')
    ax.autoscale(tight=True)    
    plt.tight_layout()
    plt.show()