# -*- coding: utf-8 -*-

import numpy                as np
import matplotlib.pyplot    as plt
import warnings

# %% constantes que ayudaran en la lectura del codigo
X, Y = 0, 1

# %%
def t2ft_H20(xnod, lado, carga, espesor):
    '''Función que convierte las fuerzas superficiales aplicadas a un elemento
    finito rectangular de 8 (serendípito) o 9 (lagrangiano) nodos a sus 
    correspondientes cargas nodales equivalentes ft    
    
    Recibe:
        xnod:  coordenadas nodales del elemento finito
          SERENDIPITO 8          LAGRANGIANO 9
          xnod = [ x1e y1e       xnod = [ x1e y1e
                  x2e y2e                x2e y2e
                   ... ...                ... ...
                  x8e y8e ]             x9e y9e ]

        lado:  arista en la que se aplica la carga, puede tomar los siguientes
               valores: 123, 345, 567, 781

        carga: fuerza distribuida en los nodos
        
               [ t1x t1y t2x t2y t3x t3y ]; % si carga se aplica sobre lado 123
               [ t3x t3y t4x t4y t5x t5y ]; % si carga se aplica sobre lado 345
               [ t5x t5y t6x t6y t7x t7y ]; % si carga se aplica sobre lado 567
               [ t7x t7y t8x t8y t1x t1y ]; % si carga se aplica sobre lado 781
    
        espesor: espesor del elemento
    '''
    
    # se definen los indices de los lados
    if   lado == 123: idx = np.array([ 1, 2, 3 ]) - 1
    elif lado == 345: idx = np.array([ 3, 4, 5 ]) - 1
    elif lado == 567: idx = np.array([ 5, 6, 7 ]) - 1
    elif lado == 781: idx = np.array([ 7, 8, 1 ]) - 1
    else: 
        raise Exception('Únicamente se permiten los lados 123, 345, 567 o 781')

    nno = xnod.shape[0]
    if nno not in (8, 9):
        raise Exception('Solo para elementos rectangulares de 8 o 9 nodos')

    # parámetros para mejorar la lectura del código
    X, Y = 0, 1
   
    # se define el número de puntos de la cuadratura y se obtienen los puntos
    # de evaluación y los pesos de la cuadratura
    n_gl       = 5
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gl)
    
    # se definen las funciones de forma unidimensionales y sus derivadas
    NN      = lambda xi: np.array([ xi*(xi-1)/2, (1+xi)*(1-xi), xi*(1+xi)/2 ])
    dNN_dxi = lambda xi: np.array([ xi - 1/2   , -2*xi        , xi + 1/2    ])
       
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

# %%
def N_H20(xi,eta,zeta):
    # Funciones de forma del EF hexaédrico serendípito isoparamétrico de 20 nodos
    return np.array([
         ((eta - 1)*(xi - 1)*(zeta - 1)*(eta + xi + zeta + 2))/8,     # N1
        -((xi**2 - 1)*(eta - 1)*(zeta - 1))/4                   ,     # N2
        -((eta - 1)*(xi + 1)*(zeta - 1)*(eta - xi + zeta + 2))/8,     # N3
         ((eta**2 - 1)*(xi + 1)*(zeta - 1))/4                   ,     # N4
        -((eta + 1)*(xi + 1)*(zeta - 1)*(eta + xi - zeta - 2))/8,     # N5
         ((xi**2 - 1)*(eta + 1)*(zeta - 1))/4                   ,     # N6
        -((eta + 1)*(xi - 1)*(zeta - 1)*(xi - eta + zeta + 2))/8,     # N7
        -((eta**2 - 1)*(xi - 1)*(zeta - 1))/4                   ,     # N8
        -((zeta**2 - 1)*(eta - 1)*(xi - 1))/4                   ,     # N9
         ((zeta**2 - 1)*(eta - 1)*(xi + 1))/4                   ,     # N10
        -((zeta**2 - 1)*(eta + 1)*(xi + 1))/4                   ,     # N11
         ((zeta**2 - 1)*(eta + 1)*(xi - 1))/4                   ,     # N12
        -((eta - 1)*(xi - 1)*(zeta + 1)*(eta + xi - zeta + 2))/8,     # N13
         ((xi**2 - 1)*(eta - 1)*(zeta + 1))/4                   ,     # N14
         ((eta - 1)*(xi + 1)*(zeta + 1)*(eta - xi - zeta + 2))/8,     # N15
        -((eta**2 - 1)*(xi + 1)*(zeta + 1))/4                   ,     # N16
         ((eta + 1)*(xi + 1)*(zeta + 1)*(eta + xi + zeta - 2))/8,     # N17
        -((xi**2 - 1)*(eta + 1)*(zeta + 1))/4                   ,     # N18
        -((eta + 1)*(xi - 1)*(zeta + 1)*(eta - xi + zeta - 2))/8,     # N19
         ((eta**2 - 1)*(xi - 1)*(zeta + 1))/4                         # N20
    ])                                                              

# %%
def dN_dxi_H20(xi,eta,zeta):
    # Derivadas de las funciones de forma con respecto a la variable "xi" del
    # EF hexaédrico serendípito isoparamétrico de 20 nodos
    return np.array([
         ((eta - 1)*(zeta - 1)*(eta + 2*xi + zeta + 1))/8,   # dN1_dxi 
        -(xi*(eta - 1)*(zeta - 1))/2                     ,   # dN2_dxi 
        -((eta - 1)*(zeta - 1)*(eta - 2*xi + zeta + 1))/8,   # dN3_dxi 
         ((eta**2 - 1)*(zeta - 1))/4                     ,   # dN4_dxi 
        -((eta + 1)*(zeta - 1)*(eta + 2*xi - zeta - 1))/8,   # dN5_dxi 
         (xi*(eta + 1)*(zeta - 1))/2                     ,   # dN6_dxi 
        -((eta + 1)*(zeta - 1)*(2*xi - eta + zeta + 1))/8,   # dN7_dxi 
        -((eta**2 - 1)*(zeta - 1))/4                     ,   # dN8_dxi 
        -((zeta**2 - 1)*(eta - 1))/4                     ,   # dN9_dxi 
         ((zeta**2 - 1)*(eta - 1))/4                     ,   # dN10_dxi
        -((zeta**2 - 1)*(eta + 1))/4                     ,   # dN11_dxi
         ((zeta**2 - 1)*(eta + 1))/4                     ,   # dN12_dxi
        -((eta - 1)*(zeta + 1)*(eta + 2*xi - zeta + 1))/8,   # dN13_dxi
         (xi*(eta - 1)*(zeta + 1))/2                     ,   # dN14_dxi
         ((eta - 1)*(zeta + 1)*(eta - 2*xi - zeta + 1))/8,   # dN15_dxi
        -((eta**2 - 1)*(zeta + 1))/4                     ,   # dN16_dxi
         ((eta + 1)*(zeta + 1)*(eta + 2*xi + zeta - 1))/8,   # dN17_dxi
        -(xi*(eta + 1)*(zeta + 1))/2                     ,   # dN18_dxi
        -((eta + 1)*(zeta + 1)*(eta - 2*xi + zeta - 1))/8,   # dN19_dxi
         ((eta**2 - 1)*(zeta + 1))/4                         # dN20_dxi
        ])

# %%
def dN_deta_H20(xi,eta,zeta):
    # Derivadas de las funciones de forma con respecto a la variable "eta" del
    # EF hexaédrico serendípito isoparamétrico de 20 nodos
    return np.array([
         ((xi - 1)*(zeta - 1)*(2*eta + xi + zeta + 1))/8,    # dN1_deta  
        -((xi**2 - 1)*(zeta - 1))/4                     ,    # dN2_deta  
        -((xi + 1)*(zeta - 1)*(2*eta - xi + zeta + 1))/8,    # dN3_deta  
         (eta*(xi + 1)*(zeta - 1))/2                    ,    # dN4_deta  
        -((xi + 1)*(zeta - 1)*(2*eta + xi - zeta - 1))/8,    # dN5_deta  
         ((xi**2 - 1)*(zeta - 1))/4                     ,    # dN6_deta  
        -((xi - 1)*(zeta - 1)*(xi - 2*eta + zeta + 1))/8,    # dN7_deta  
        -(eta*(xi - 1)*(zeta - 1))/2                    ,    # dN8_deta  
        -((zeta**2 - 1)*(xi - 1))/4                     ,    # dN9_deta  
         ((zeta**2 - 1)*(xi + 1))/4                     ,    # dN10_deta 
        -((zeta**2 - 1)*(xi + 1))/4                     ,    # dN11_deta 
         ((zeta**2 - 1)*(xi - 1))/4                     ,    # dN12_deta 
        -((xi - 1)*(zeta + 1)*(2*eta + xi - zeta + 1))/8,    # dN13_deta 
         ((xi**2 - 1)*(zeta + 1))/4                     ,    # dN14_deta 
         ((xi + 1)*(zeta + 1)*(2*eta - xi - zeta + 1))/8,    # dN15_deta 
        -(eta*(xi + 1)*(zeta + 1))/2                    ,    # dN16_deta 
         ((xi + 1)*(zeta + 1)*(2*eta + xi + zeta - 1))/8,    # dN17_deta 
        -((xi**2 - 1)*(zeta + 1))/4                     ,    # dN18_deta 
        -((xi - 1)*(zeta + 1)*(2*eta - xi + zeta - 1))/8,    # dN19_deta 
         (eta*(xi - 1)*(zeta + 1))/2                         # dN20_deta 
    ])

# %%
def dN_dzeta_H20(xi,eta,zeta):
    # Derivadas de las funciones de forma con respecto a la variable "zeta" del
    # EF hexaédrico serendípito isoparamétrico de 20 nodos
    return np.array([
          ((eta - 1)*(xi - 1)*(eta + xi + 2*zeta + 1))/8,    # dN1_dzeta 
         -((xi**2 - 1)*(eta - 1))/4                     ,    # dN2_dzeta 
         -((eta - 1)*(xi + 1)*(eta - xi + 2*zeta + 1))/8,    # dN3_dzeta 
          ((eta**2 - 1)*(xi + 1))/4                     ,    # dN4_dzeta 
         -((eta + 1)*(xi + 1)*(eta + xi - 2*zeta - 1))/8,    # dN5_dzeta 
          ((xi**2 - 1)*(eta + 1))/4                     ,    # dN6_dzeta 
         -((eta + 1)*(xi - 1)*(xi - eta + 2*zeta + 1))/8,    # dN7_dzeta 
         -((eta**2 - 1)*(xi - 1))/4                     ,    # dN8_dzeta 
         -(zeta*(eta - 1)*(xi - 1))/2                   ,    # dN9_dzeta 
          (zeta*(eta - 1)*(xi + 1))/2                   ,    # dN10_dzeta
         -(zeta*(eta + 1)*(xi + 1))/2                   ,    # dN11_dzeta
          (zeta*(eta + 1)*(xi - 1))/2                   ,    # dN12_dzeta
         -((eta - 1)*(xi - 1)*(eta + xi - 2*zeta + 1))/8,    # dN13_dzeta
          ((xi**2 - 1)*(eta - 1))/4                     ,    # dN14_dzeta
          ((eta - 1)*(xi + 1)*(eta - xi - 2*zeta + 1))/8,    # dN15_dzeta
         -((eta**2 - 1)*(xi + 1))/4                     ,    # dN16_dzeta
          ((eta + 1)*(xi + 1)*(eta + xi + 2*zeta - 1))/8,    # dN17_dzeta
         -((xi**2 - 1)*(eta + 1))/4                     ,    # dN18_dzeta
         -((eta + 1)*(xi - 1)*(eta - xi + 2*zeta - 1))/8,    # dN19_dzeta
          ((eta**2 - 1)*(xi - 1))/4                          # dN20_dzeta
        ])

def matriz_extrapolacion_esfuerzos_H20():
    return np.array([
       [(3*3**(1/2))/4 + 5/4,   - 3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4],
       [    3**(1/2)/4 + 1/2,                 -1/4,                 -1/4,     1/2 - 3**(1/2)/4,     3**(1/2)/4 + 1/2,                 -1/4,                 -1/4,     1/2 - 3**(1/2)/4],
       [  - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4, (3*3**(1/2))/4 + 5/4,   - 3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4],
       [                -1/4,     1/2 - 3**(1/2)/4,                 -1/4,     1/2 - 3**(1/2)/4,     3**(1/2)/4 + 1/2,                 -1/4,     3**(1/2)/4 + 1/2,                 -1/4],
       [    3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4,   - 3**(1/2)/4 - 1/4],
       [                -1/4,     1/2 - 3**(1/2)/4,     3**(1/2)/4 + 1/2,                 -1/4,                 -1/4,     1/2 - 3**(1/2)/4,     3**(1/2)/4 + 1/2,                 -1/4],
       [  - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4],
       [    3**(1/2)/4 + 1/2,                 -1/4,     3**(1/2)/4 + 1/2,                 -1/4,                 -1/4,     1/2 - 3**(1/2)/4,                 -1/4,     1/2 - 3**(1/2)/4],
       [    3**(1/2)/4 + 1/2,     3**(1/2)/4 + 1/2,                 -1/4,                 -1/4,                 -1/4,                 -1/4,     1/2 - 3**(1/2)/4,     1/2 - 3**(1/2)/4],
       [                -1/4,                 -1/4,     1/2 - 3**(1/2)/4,     1/2 - 3**(1/2)/4,     3**(1/2)/4 + 1/2,     3**(1/2)/4 + 1/2,                 -1/4,                 -1/4],
       [    1/2 - 3**(1/2)/4,     1/2 - 3**(1/2)/4,                 -1/4,                 -1/4,                 -1/4,                 -1/4,     3**(1/2)/4 + 1/2,     3**(1/2)/4 + 1/2],
       [                -1/4,                 -1/4,     3**(1/2)/4 + 1/2,     3**(1/2)/4 + 1/2,     1/2 - 3**(1/2)/4,     1/2 - 3**(1/2)/4,                 -1/4,                 -1/4],
       [  - 3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4,     3**(1/2)/4 - 1/4],
       [                -1/4,     3**(1/2)/4 + 1/2,     1/2 - 3**(1/2)/4,                 -1/4,                 -1/4,     3**(1/2)/4 + 1/2,     1/2 - 3**(1/2)/4,                 -1/4],
       [    3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4, 5/4 - (3*3**(1/2))/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4],
       [    1/2 - 3**(1/2)/4,                 -1/4,     1/2 - 3**(1/2)/4,                 -1/4,                 -1/4,     3**(1/2)/4 + 1/2,                 -1/4,     3**(1/2)/4 + 1/2],
       [5/4 - (3*3**(1/2))/4,     3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4],
       [    1/2 - 3**(1/2)/4,                 -1/4,                 -1/4,     3**(1/2)/4 + 1/2,     1/2 - 3**(1/2)/4,                 -1/4,                 -1/4,     3**(1/2)/4 + 1/2],
       [    3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4, (3*3**(1/2))/4 + 5/4, 5/4 - (3*3**(1/2))/4,     3**(1/2)/4 - 1/4,     3**(1/2)/4 - 1/4,   - 3**(1/2)/4 - 1/4],
       [                -1/4,     3**(1/2)/4 + 1/2,                 -1/4,     3**(1/2)/4 + 1/2,     1/2 - 3**(1/2)/4,                 -1/4,     1/2 - 3**(1/2)/4,                 -1/4]])
                       
def gausslegendre_quad_hexa(n_gl):
    x, w = np.polynomial.legendre.leggauss(n_gl)

    x_gl = np.zeros((n_gl**3,3))
    w_gl = np.zeros(n_gl**3)

    r = -1
    for i in range(n_gl):
        for j in range(n_gl):
            for k in range(n_gl):
                r += 1
                x_gl[r,:] = np.array([ x[i], x[j], x[k] ])
                w_gl[r]   = w[i]*w[j]*w[k]

    return x_gl, w_gl


'''
# %% Se dibuja la malla de elementos finitos
cg = np.zeros((nef,2))  # almacena el centro de gravedad de los EF
plt.figure()
for e in range(nef):
    # se dibujan las aristas
    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
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
'''

'''
# %% Dibujo la malla de elementos finitos y las deformada de esta
delta  = np.reshape(a, (nno,2))
escala = 50000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

plt.figure()
for e in range(nef):
   nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
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
'''