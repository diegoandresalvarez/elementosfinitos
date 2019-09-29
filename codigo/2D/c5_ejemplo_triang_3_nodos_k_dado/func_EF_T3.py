# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# se definen algunas constantes
X, Y = 0, 1
NL1, NL2, NL3 = 0, 1, 2

def t2ft_T3(xnod, lado, carga, espesor):
    '''Esta función convierte las fuerzas superficiales aplicadas a un elemento
    finito triangular de 3 nodos a sus correspondientes cargas nodales
    equivalentes ft

    xnod = np.array([[x1e, y1e,
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
        L12 = np.hypot(xnod[NL1,X] - xnod[NL2,X], xnod[NL1,Y] - xnod[NL2,Y])

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
        L23 = np.hypot(xnod[NL2,X] - xnod[NL3,X], xnod[NL2,Y] - xnod[NL3,Y])

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
        L31 = np.hypot(xnod[NL3,X] - xnod[NL1,X], xnod[NL3,Y] - xnod[NL1,Y])

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

