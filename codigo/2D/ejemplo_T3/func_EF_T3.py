# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors

# se definen algunas constantes
X, Y = 0, 1
NL1, NL2, NL3 = 0, 1, 2

# variables globales que se heredarán del programa principal
xnod = None
LaG = None
cg = None
interpolar_colores = None

# variable global que calculará compartir_variables()
num_elem_ady = None

def compartir_variables(xnod_, LaG_, cg_, interpolar=False):
    '''Importa variables globales del programa principal a este módulo
    '''
    global xnod, LaG, cg, interpolar_colores, num_elem_ady
    xnod, LaG, cg = xnod_, LaG_, cg_
    interpolar_colores = interpolar
    if interpolar_colores:
        nno = xnod.shape[0]
        nef = LaG.shape[0]

        # El array "num_elem_ady" contabiliza el número de EFs adyacentes a un
        # nodo dado, con el objeto de hacer luego el alisado de los esfuerzos y
        # las deformaciones
        num_elem_ady = np.zeros(nno)
        for e in range(nef):
            num_elem_ady[LaG[e,:]] += 1

def t2ft_T3(xnod_EF, lado, carga, espesor):
    '''Convierte las fuerzas superficiales aplicadas a un EF triangular de 3
    nodos a sus correspondientes cargas nodales equivalentes ft

    xnod_EF = np.array([[x1e, y1e],
                        [x2e, y2e],
                        [x3e, y3e]])

    lado = 12, 23, 31

    carga = [ t1x, t1y, t2x, t2y ]   # si la carga se aplica sobre el lado 12
            [ t2x, t2y, t3x, t3y ]   # si la carga se aplica sobre el lado 23
            [ t3x, t3y, t1x, t1y ]   # si la carga se aplica sobre el lado 31
    '''

    if lado == 12:
        # Fuerzas sobre el lado 12
        # Se calcula la longitud del lado 12
        L12 = np.hypot(xnod_EF[NL1,X] - xnod_EF[NL2,X], xnod_EF[NL1,Y] - xnod_EF[NL2,Y])

        # Fuerzas distribuidas aplicadas en los nodos 1 y 2 locales
        t1x, t1y, t2x, t2y = carga
        ft = espesor*np.array([ (L12*(2*t1x + t2x))/6,
                                (L12*(2*t1y + t2y))/6,
                                (L12*(t1x + 2*t2x))/6,
                                (L12*(t1y + 2*t2y))/6,
                                0,
                                0                      ])
    elif lado == 23:
        # Fuerzas sobre el lado 23
        # Se calcula la longitud del lado 23
        L23 = np.hypot(xnod_EF[NL2,X] - xnod_EF[NL3,X], xnod_EF[NL2,Y] - xnod_EF[NL3,Y])

        # Fuerzas distribuidas aplicadas en los nodos 2 y 3 locales
        t2x, t2y, t3x, t3y = carga
        ft = espesor*np.array([ 0,
                                0,
                                (L23*(2*t2x + t3x))/6,
                                (L23*(2*t2y + t3y))/6,
                                (L23*(t2x + 2*t3x))/6,
                                (L23*(t2y + 2*t3y))/6 ])
    elif lado == 31:
        # Fuerzas sobre el lado 31
        # Se calcula la longitud del lado 31
        L31 = np.hypot(xnod_EF[NL3,X] - xnod_EF[NL1,X], xnod_EF[NL3,Y] - xnod_EF[NL1,Y])

        # Fuerzas distribuidas aplicadas en los nodos 3 y 1 locales
        t3x, t3y, t1x, t1y = carga
        ft = espesor*np.array([ (L31*(2*t1x + t3x))/6,
                                (L31*(2*t1y + t3y))/6,
                                0,
                                0,
                                (L31*(t1x + 2*t3x))/6,
                                (L31*(t1y + 2*t3y))/6 ])
    else:
        raise Exception('Unicamente se permiten los lados 12, 23 o 31')

    return ft

'''
# NOTA: Los vectores se calcularon con el siguiente programa de MATLAB:
# Elemento triangular de 3 nodos

clear
clc
syms N1 N2 N3
syms t1x t1y t2x t2y t3x t3y
syms xi L12 L23 L31

N = [ N1 0   N2 0   N3 0
      0  N1  0  N2  0  N3 ];

t = [ t1x t1y t2x t2y t3x t3y ].';

res = N.'*N*t;

res12 = res;
res12 = subs(res12, N1, (1-xi)/2);
res12 = subs(res12, N2, (1+xi)/2);
res12 = subs(res12, N3, 0);
ds_dxi = L12/2;
f12 = int(res12*ds_dxi, xi, -1 ,1)

res23 = res;
res23 = subs(res23, N2, (1-xi)/2);
res23 = subs(res23, N3, (1+xi)/2);
res23 = subs(res23, N1, 0);
ds_dxi = L23/2;
f23 = int(res23*ds_dxi, xi, -1 ,1)

res31 = res;
res31 = subs(res31, N3, (1-xi)/2);
res31 = subs(res31, N1, (1+xi)/2);
res31 = subs(res31, N2, 0);
ds_dxi = L31/2;
f31 = int(res31*ds_dxi, xi, -1 ,1)
'''

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import warnings

def plot_esf_def(variable, titulo, angulo=None):
    '''Grafica los esfuerzos y las deformaciones para la malla de EFs dada.

    Es muy importante tener en cuenta que previamente se debió haber ejecutado
    la función "compartir_variables()" con el objeto de incluir en el presente
    módulo algunas variables globales como xnod y LaG que serán utilizas en este
    programa.

    Uso:
        plot_esf_def(variable, titulo, angulo=None):

    variable: es la variable que se quiere graficar
    titulo:   del gráfico
    angulo:   para los esfuerzos principales s1 y s2 y tmax
    '''

    # se determina el número de nodos y de EFs
    nno = xnod.shape[0]
    nef = LaG.shape[0]

    # se inicializa el lienzo
    fig, ax = plt.subplots()

    # y se hace el gráfico respectivo

    if interpolar_colores:
        # se promedian los esfuerzos y las deformaciones en los nodos de modo
        # que en la variable "var" se encuentren los valores alisados del
        # esfuerzo o de la deformación a graficar
        var = np.zeros(nno)
        for e in range(nef):
            var[LaG[e,:]] += variable[e]
        var /= num_elem_ady

        # se encuentra el máximo en valor absoluto para ajustar el colorbar()
        val_max = np.max(np.abs(var))
        #val_max = np.max(var)
        #val_min = np.min([np.min(var), -np.finfo(float).eps])

        # se grafica la malla de EFS, los colores en cada triángulo y las curvas
        # de nivel
        ax.triplot        (xnod[:,X], xnod[:,Y], LaG, lw=0.5, color='gray')
        # norm = mcolors.DivergingNorm(vmin=val_min, vmax = val_max, vcenter=0)
        im = ax.tripcolor (xnod[:,X], xnod[:,Y], LaG, var, cmap='bwr',
                                shading='gouraud', vmin=-val_max, vmax=val_max) #norm=norm)
        ax.tricontour(xnod[:,X], xnod[:,Y], LaG, var, 20)

        # a veces sale un warning, simplemente porque no existe la curva 0
        warnings.filterwarnings("ignore")
        ax.tricontour(xnod[:,X], xnod[:,Y], LaG, var, levels=[0], linewidths=3)
        warnings.filterwarnings("default")

        fig.colorbar(im, ax=ax, format='%6.3g')
    else:
        # cada EF se especifica como un polígono
        patches = [ Polygon(np.c_[xnod[LaG[e,:],X], xnod[LaG[e,:],Y]], closed=True)
                                                            for e in range(nef) ]

        # se encuentra el máximo en valor absoluto para ajustar el colorbar()
        val_max = np.max(np.abs(variable))

        # se crean cada uno de los parches, con su valor asociado y colormap respectivo
        p = PatchCollection(patches, array=variable, cmap='bwr')

        ax.add_collection(p)

        p.set_clim(vmin=-val_max, vmax=val_max)
        fig.colorbar(p, ax=ax, format='%6.3g')

    # se grafican las líneas que indiquen direcciones de los esfuerzos
    if angulo is not None:
       for ang in angulo:
            ax.quiver(cg[:,X], cg[:,Y],
                variable*np.cos(ang), variable*np.sin(ang),
                headwidth=0, headlength=0, headaxislength=0,
                pivot='middle')

    # se especifican los ejes y el título, y se colocan los ejes iguales
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(titulo, fontsize=20)

    # esto es como un "axis tight" de matlab
    ax.set_aspect('equal')
    ax.autoscale(tight=True)

    plt.show()

'''
PENDIENTE:
* hacer que al mover el mouse, se muestren en el título los esfuerzos de triángulo
'''
