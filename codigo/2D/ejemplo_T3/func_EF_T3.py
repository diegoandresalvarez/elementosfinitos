# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# se definen algunas constantes
X, Y = 0, 1
NL1, NL2, NL3 = 0, 1, 2
xnod = None
LaG = None
cg = None
interpolar_colores = None
num_elem_ady = None

def compartir_variables(xnod_, LaG_, cg_, interpolar=False):
    global xnod
    global LaG
    global cg
    global interpolar_colores
    global num_elem_ady
    xnod = xnod_
    LaG  = LaG_
    cg   = cg_
    interpolar_colores = interpolar
    if interpolar_colores:
        nno = xnod.shape[0]
        nef = LaG.shape[0]
        num_elem_ady = np.zeros(nno)  # numero de elementos adyacentes
        for e in range(nef):
            num_elem_ady[LaG[e,:]] += 1

def t2ft_T3(xnod_EF, lado, carga, espesor):
    '''Esta función convierte las fuerzas superficiales aplicadas a un elemento
    finito triangular de 3 nodos a sus correspondientes cargas nodales
    equivalentes ft

    xnod_EF = np.array([[x1e, y1e,
                     x2e, y2e,
                     x3e, y3e]]
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
dx_dxi = L12/2;
f12 = int(res12*dx_dxi, xi, -1 ,1)

res23 = res;
res23 = subs(res23, N2, (1-xi)/2);
res23 = subs(res23, N3, (1+xi)/2);
res23 = subs(res23, N1, 0);
dx_dxi = L23/2;
f23 = int(res23*dx_dxi, xi, -1 ,1)

res31 = res;
res31 = subs(res31, N3, (1-xi)/2);
res31 = subs(res31, N1, (1+xi)/2);
res31 = subs(res31, N2, 0);
dx_dxi = L31/2;
f31 = int(res31*dx_dxi, xi, -1 ,1)
'''

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
#from matplotlib import colors


def plot_esf_def(variable, titulo, angulo=None):
    '''FALTA
    '''
    # se determina el número de elementos finitos
    nno = xnod.shape[0]
    nef = LaG.shape[0]    
    
    if interpolar_colores:
        # se promedian los esfuerzos y las deformaciones en los nodos
        variable_ = np.zeros(nno)
        for e in range(nef):
            variable_[LaG[e,:]] += variable[e]
        variable_ /= num_elem_ady
        fig, ax = plt.subplots()
        ax.triplot(xnod[:,X], xnod[:,Y], LaG, lw=0.5, color='white')
        ax.tripcolor(xnod[:,X], xnod[:,Y], LaG, variable_, cmap=cm.jet, shading='gouraud')
        ax.tricontour(xnod[:,X], xnod[:,Y], LaG, variable_, 20, 
            colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
            linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])                
        #plt.colorbar() #ax=ax,format='%.3g')        
    else:
        fig, ax = plt.subplots()
    
        # cada EF se especifica como un polígono
        patches = [ Polygon(np.c_[xnod[LaG[e,:],X], xnod[LaG[e,:],Y]], closed=True) 
                                                            for e in range(nef) ]

        #divnorm = colors.DivergingNorm(vmin=min(variable), vcenter=0., vmax=max(variable))                                                           
        #p = PatchCollection(patches, cmap=cm.bwr, norm=divnorm) #, alpha=0.4)
        p = PatchCollection(patches, cmap=cm.jet) #, alpha=0.4)
    
        # y se les especifica el color asociado
        p.set_array(variable)

        ax.add_collection(p)
        # se crea la barra de colores con la misma altura del gráfico, a la derecha 
        # del mismo ( a la derecha de ax). El ancho será el 3% del ancho de ax y el 
        # espacio entre cax y ax será de 0.3 pulgadas
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.3)
        plt.colorbar(p, cax=cax, format='%.3g')        
    ax.autoscale(enable=True, tight=True)

    # se especifican los ejes y el título, y se colocan los ejes iguales
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title(titulo)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout()

    # se grafican las líneas que indiquen direcciones de los esfuerzos    
    if angulo is not None:
       esc = 1.1 # escala para graficar las flechas
       for ang in angulo:
            ax.quiver(cg[:,X], cg[:,Y],
                variable*np.cos(ang), variable*np.sin(ang),
                headwidth=0)
            ax.quiver(cg[:,X], cg[:,Y],
                variable*np.cos(ang+np.pi), variable*np.sin(ang+np.pi),
                headwidth=0)         
    plt.show()        
 
'''
PENDIENTE:
* hacer que al mover el mouse, se muestren en el título los esfuerzos de triángulo
* comentar mejor este archivo 
* mejorar el README.md
'''
