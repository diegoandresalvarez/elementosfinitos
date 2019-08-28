# -*- coding: utf-8 -*-

from math import atan2
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

def c1_dibujar_barra_deformada_portico(A, E, I, x1,y1, x2,y2, qxloc,qyloc, 
                                       qe, ae, esc_def, esc_faxial, esc_V, esc_M):
    '''
    Esta función dibuja el elemento de pórtico deformado junto con sus 
    respectivos diagramas de fuerza axial, fuerza cortante y momento flector.

    El diagrama de momento flector se grafica en el lado opuesto de la fibra
    a tracción

    PARAMETROS DE ENTRADA (junto con algunos ejemplos):
    A = area
    E = E
    I = Ix local
    (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
    qxloc = lambda x : x**2 # carga en la dir. del eje x local (function handle)
    qyloc = lambda x : 0    # carga en la dir. del eje y local (function handle)
    qe = np.array(
         [ 0.01,         # U1, V1, M1 reacciones del nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,         # U2, V2, M2 reacciones del nodo 2 en coord. locales
           0.02,
          -0.07 ])
    ae = np.array(
         [ 0.01,         # u1, v1, t1 desplazamientos nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,         # u2, v2, t2 desplazamientos nodo 2 en coord. locales
           0.02,
          -0.07 ])
    esc_def    = 10     # escalamiento de la deformada
    esc_faxial = 10     # escalamiento del diagrama de axiales
    esc_V      = 10     # escalamiento del diagrama de cortantes
    esc_M      = 10     # escalamiento del diagrama de momentos
    '''

    #%% se definen algunas constantes
    X = 0; Y = 1; X1 = 0; Y1 = 1; M1 = 2; X2 = 3; Y2 = 4; M2 = 5
    
    #%% se define la ecuación diferencial asociada
    def ecuacion_diferencial(x,y):
        # aquí se implementa la ecuación diferencial para vigas de material
        # homogéneo y sección transversal constante (A, E, I, qx, qy las 
        # provee la función exterior)
        #      d^4 v(x)
        # E I ---------- = q(x)
        #        dx^4
        #
        #      d^2 u(x)
        # A E ---------- = -b(x)
        #        dx^2
        
        m = x.shape[0]
        dydx = np.zeros((6,m))
        #           y[0,:]         = v
        dydx[0,:] = y[1,:]       # = theta
        dydx[1,:] = y[2,:]/(E*I) # = M/(EI)
        dydx[2,:] = y[3,:]       # = V
        dydx[3,:] = qyloc(x)     # = qyloc
        dydx[4,:] = y[5,:]/(A*E) # = u
        dydx[5,:] = -qxloc(x)    # = faxial
          
        return dydx
    
    
    #%% ------------------------------------------------------------------------
    
    def condiciones_de_apoyo(YL,YR):
        # condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
        u1  = 0; v1 = 1; t1 = 2; u2 = 3; v2 = 4; t2   = 5
        v_  = 0; t_ = 1;                 u_ = 4
        #                M_ = 2; V_ = 3;        fax_ = 5
        return np.array([ 
                    # YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
                    YL[u_] - ae[u1],        # uloc(0)     = u1
                    YL[v_] - ae[v1],        # vloc(0)     = v1
                    YL[t_] - ae[t1],        # thetaloc(0) = t1
                    YR[u_] - ae[u2],        # uloc(L)     = u2
                    YR[v_] - ae[v2],        # vloc(L)     = v2
                    YR[t_] - ae[t2] ])      # thetaloc(L) = t2
        
    #%% resolver la ecuación diferencial
    npuntos = 1001
    xinit = np.linspace(0, np.hypot(x2-x1, y2-y1), npuntos)
    yinit = np.zeros((6,npuntos))
    sol   = solve_bvp(ecuacion_diferencial, condiciones_de_apoyo, xinit, yinit)
    #https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.integrate.solve_bvp.html#scipy.integrate.solve_bvp
    
    #%% Calculos intermedios
    s     = sol.x
    axial = sol.y[5,:]            # Fuerza axial [kN]
    V     = sol.y[3,:]            # Fuerza cortante [kN]
    M     = sol.y[2,:]            # Momento flector [kN/m]
    u     = sol.y[4,:]            # Desplazamiento horizontal de la viga [m]
    v     = sol.y[0,:]            # Desplazamiento vertical de la viga [m]
    #theta = np.arctan(sol.y[1,:]) # Angulo de giro  [rad]
    
    # rotación de la solución antes de dibujar
    ang = atan2(y2-y1, x2-x1)
    T   = np.array([[np.cos(ang),  -np.sin(ang)], # matriz de rotación
                    [np.sin(ang),   np.cos(ang)]])
      
    #%% Dibujar de deformada
    plt.figure(2)
    pos = T @ np.array([s + esc_def*u,
                            esc_def*v])
    
    xx = pos[X,:] + x1
    yy = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-', xx, yy, 'r-', linewidth=2)
    
    #%% Dibujar los diagramas de fuerza axial 
    plt.figure(3)
    pos = T @ np.array([s               ,
                        esc_faxial*axial]) # escalamiento del diagrama
    
    ss = pos[X,:] + x1
    aa = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-', np.r_[x1, ss, x2], np.r_[y1, aa, y2], 'r-', linewidth=2)
    plt.text(ss[ 0], aa[ 0], str(-qe[X1]))
    plt.text(ss[-1], aa[-1], str(+qe[X2]))
    
    #%% Dibujar los diagramas de fuerza cortante
    plt.figure(4)
    pos = T @ np.array([s      ,
                        esc_V*V]) # escalamiento del diagrama
    
    ss = pos[X,:] + x1
    vv = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-', np.r_[x1, ss, x2], np.r_[y1, vv, y2], 'r-', linewidth=2)
    plt.text(ss[ 0], vv[ 0], str(+qe[Y1]))
    plt.text(ss[-1], vv[-1], str(-qe[Y2]))
    
    #%% Dibujar los diagramas de momento flector
    plt.figure(5)
    pos = T @ np.array([s      ,
                        esc_M*M]) # escalamiento del diagrama    
    
    ss = pos[X,:] + x1
    mm = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-', np.r_[x1, ss, x2], np.r_[y1, mm, y2], 'r-', linewidth=2)
    plt.text(ss[ 0], mm[ 0], str(-qe[M1]))
    plt.text(ss[-1], mm[-1], str(+qe[M2]))

#%% bye, bye!
