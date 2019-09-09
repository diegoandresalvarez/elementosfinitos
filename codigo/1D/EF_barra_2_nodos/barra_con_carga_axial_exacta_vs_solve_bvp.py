# -*- coding: utf-8 -*-

# %% definición del problema
# Calcule los desplazamientos y las reacciones en el empotramiento
# de la barra mostrada resolviendo la ecuación diferencial numéricamente con
# la funcion bvp4c
#
# | b (carga distribuida de magnitud b)
# |->->->->->->->->->->->->->->->->
# |================================o--> P (carga puntual P en extremo derecho)
# |<----longitud L de la barra---->|    el area transversal de la barra es A

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
# %matplotlib --list
# %matplotlib auto
# %matplotlib inline

# %% defino las variables
E   = 200e9     # Pa          # módulo de elasticidad de la barra
A   = (0.01)**2 # m^2         # área transversal de la barra
L   = 2         # m           # longitud de la barra
b   = 1000      # N/m         # fuerza axial aplicada sobre cada EF
P   = 250       # N           # carga nodal al final de la barra

## Solución de la ecuacion diferencial

# Solución numérica usando bvp4c (boundary value problem - MATLAB)
#   d /           du(x)  \
# ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0)=0
#  dx \            dx    /                                faxial(L) = P
#
# bvp4c is a finite difference code that implements the three-stage Lobatto
# IIIa formula. This is a collocation formula and the collocation polynomial
# provides a C1-continuous solution that is fourth-order accurate uniformly 
# in [a,b]. Mesh selection and error control are based on the residual of 
# the continuous solution. 

# En el caso mas general E, A y b son funciones. Escriba aquí las funciones
# como tal en caso de tener un caso mas general
EE = lambda x : E
AA = lambda x : A
bb = lambda x : b

# Por favor, antes de continuar mirar la ayuda de MATLAB de los comandos
# bvpinit, bvp4c, deval

# %% se define la ecuación diferencial asociada
def ecuacion_diferencial(x,y):
    # aquí se implementa la ecuación diferencial para vigas de material
    # homogéneo y sección transversal constante (A, E, I, qx, qy las 
    # provee la función exterior)

    # Se define la ecuacion diferencial, expresada como un sistema de dos
    # ecuaciones diferenciales
    # y(1) = u(x)
    # y(2) = faxial(x)    

    #      d^2 u(x)
    # A E ---------- = -b(x)
    #        dx^2
    
    m = x.shape[0]
    dydx = np.zeros((2,m))
    dydx[0,:] = y[1,:]/(EE(x)*AA(x)) # = u
    dydx[1,:] = -bb(x)               # = faxial
        
    return dydx

# Se definen las condiciones de frontera
# ya = condiciones de frontera del lado izquierdo (x=0)
#     ya(1) = u(x=0)          ya(2) = faxial(x=0)
# yb = condiciones de frontera del lado derecho   (x=L)
#     yb(1) = u(x=L)          yb(2) = faxial(x=L)
cond_frontera = @ (ya,yb) [ ya(1)         # u(x=0)      = 0 (desplazamiento)
                            yb(2) - P ];  # faxial(x=L) = P (carga axial)

# Solución tentativa de la ecuacion diferencial
x = np.linspace(0,L,30)         # 30 puntos uniformemente distrib. entre 0 y L
sol_inicial = bvpinit(x,[0 0]); # el [ 0 0 ] hace sol_inicial.y = zeros(2,30)

# Solución como tal de la ecuacion diferencial
sol = bvp4c(sist_eq_dif, cond_frontera, sol_inicial);

# Evaluar la respuesta en los puntos x
y = deval(sol,x)

## Solución analítica
u_exacta      = lambda x : (-b*x.^2/2 + (P + b*L)*x)/(E*A); # desplazamiento
faxial_exacta = lambda x : (P + b*(L-x));                   # carga axial

## Grafico la solución analítica y la solución por el la funcion bvp4c
# 1) grafico los desplazamientos de la barra
plt.figure                        # cree un nuevo lienzo
plt.subplot(2,1,1)                # grafique en la parte superior (1) del lienzo
xx = np.linspace(0, L, 100)       # 100 puntos equidistantes entre 0 y L
plt.plot(xx, u_exacta(xx), 'r', label='solución exacta de $u(x)$')
plt.plot(x, y(1,:), 'bx',       label='solución aproximada por bvp4c()')
plt.title('Comparación de la solución analítica vs la función bvp4c para el desplazamiento')
plt.xlabel('Eje X (m)')           # titulo del eje X
plt.ylabel('Desplazamiento (m)')  # titulo del eje Y
plt.legend(loc='lower right')

error_en_u = 100*max(abs((u_exacta(x) - y(1,:))./u_exacta(x)));
print('Maximo porcentaje de error en el calculo del desplazamiento = %g%%\n' % error_en_u)

# 2) grafico la carga axial de la barra
faxial_exacta = lambda x : P + b*(L-x) # solución analítica para la carga axial

plt.subplot(2,1,2)                # grafique en la parte inferior (2) del lienzo
plt.plot(xx, faxial_exacta(xx), 'r', label='solución exacta de $f_{axial}(x)$')
for e in np.arange(nef):          # ciclo sobre todas los elementos finitos
    plt.plot([xnod[e], xnod[e+1]], [faxial[e], faxial[e]], 'b.-',
                                  label='solución por el MEF' if e == 0 else "")

plt.title('Comparación de la solución analítica con el MEF para la carga axial')
plt.xlabel('Eje X (m)')           # título del eje X
plt.ylabel('Carga axial (N)')     # título del eje Y
plt.legend(loc='upper right')
plt.show()

#bye, bye!!!









# 2) grafico la carga axial de la barra
subplot(2,1,2);                    # grafique en la parte inferior (2) del lienzo
plot(x, faxial_exacta(x), 'r');    # grafico solución analítica
hold on;                           # no borre el lienzo
plot(x, y(2,:), 'bx');             # grafico solución por bvp4c
title('Comparacion de la solucion analitica vs la funcion bvp4c para la carga axial');
xlabel('Eje X (m)')                # título del eje X
ylabel('Carga axial (N)')          # título del eje Y
legend('solucion exacta','solucion por bvp4c', 'Location','NorthEast');

error_en_faxial = 100*max(abs((faxial_exacta(x) - y(2,:))./faxial_exacta(x)));
fprintf('Maximo porcentaje de error en el calculo de la fuerza axial = #g##\n', error_en_faxial);

## bye, bye!!!







    
    

    # %% y sus respectivas condiciones de apoyo    
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
    
