# -*- coding: utf-8 -*-

# %% Definición del problema
# Calcule los desplazamientos y las reacciones en el empotramiento
# de la barra mostrada resolviendo la ecuación diferencial numéricamente con
# la funcion solve_bvp
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

# %% Defino las variables
E   = 200e9     # Pa          # módulo de elasticidad de la barra
A   = (0.01)**2 # m^2         # área transversal de la barra
L   = 2         # m           # longitud de la barra
b   = 1000      # N/m         # fuerza axial aplicada sobre cada EF
P   = 250       # N           # carga nodal al final de la barra

# %% Solución de la ecuacion diferencial

# Solución numérica usando solve_bvp (boundary value problem)
#   d /           du(x)  \
# ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0)=0
#  dx \            dx    /                                faxial(L) = P
#

# En el caso mas general E, A y b son funciones. Escriba aquí las funciones
# como tal en caso de tener un caso mas general. 
# El vector resultante debe tener el mismo tamaño que x
EE = lambda x : np.tile(E, x.shape)
AA = lambda x : np.tile(A, x.shape)
bb = lambda x : np.tile(b, x.shape)

# Se define la ecuación diferencial expresada como un sistema de dos 
# ecuaciones diferenciales
# y[0] = u(x)
# y[1] = faxial(x)    
sist_eq_dif = lambda x,y : np.vstack([ y[1,:]/(EE(x)*AA(x)),    # = u
                                      -bb(x)                 ]) # = faxial

# Se definen las condiciones de frontera
# y_izq = condiciones de frontera del lado izquierdo (x=0)
# y_izq[0] = u(x=0)          y_izq[1] = faxial(x=0)
# y_der = condiciones de frontera del lado derecho   (x=L)
# y_der[0] = u(x=L)          y_der[1] = faxial(x=L)

# Se definen las siguientes constantes para facilitar la lectura del código
u, faxial = 0, 1

cond_frontera = lambda y_izq,y_der : \
                        [ y_izq[u],           # u(x=0)     = 0 (desplazamiento)
                          y_der[faxial] - P ] # faxial(x=L) = P (carga axial)    

# Solución tentativa de la ecuacion diferencial
x         = np.linspace(0, L, 30)   # 30 puntos uniformemente distrib. entre 0 y L
y_inicial = np.zeros((2, x.size))

# Solución como tal de la ecuacion diferencial
sol = solve_bvp(sist_eq_dif, cond_frontera, x, y_inicial)

# Evaluar la respuesta en los puntos x
y = sol.y

# %% Solución analítica
u_exacta      = lambda x : (-b*x**2/2 + (P + b*L)*x)/(E*A); # desplazamiento
faxial_exacta = lambda x : (P + b*(L-x));                   # carga axial

# %% Se reportan los errores en el calculo
error_en_u = max(abs(u_exacta(x) - y[u,:]))
print(f'Maximo error en el calculo del desplazamiento = {error_en_u} m')

error_en_faxial = max(abs(faxial_exacta(x) - y[faxial,:]))
print(f'Maximo error en el calculo de la fuerza axial = {error_en_faxial} N')

# %% Grafico la solución analítica y la solución por el la funcion solve_bvp
plt.figure                        # cree un nuevo lienzo

# 1) grafico los desplazamientos de la barra
plt.subplot(2,1,1)                # grafique en la parte superior (1) del lienzo
plt.plot(x, u_exacta(x), 'r', label='solución exacta de $u(x)$')
plt.plot(x, y[u,:], 'bx',     label='solución aproximada por solve_bvp()')
plt.title('Comparación de la solución analítica vs la función solve_bvp() para el desplazamiento')
plt.xlabel('Eje X (m)')           # titulo del eje X
plt.ylabel('Desplazamiento (m)')  # titulo del eje Y
plt.legend(loc='lower right')

# 2) grafico la carga axial de la barra
plt.subplot(2,1,2)                # grafique en la parte inferior (2) del lienzo
plt.plot(x, faxial_exacta(x), 'r', label='solución exacta de $f_{axial}(x)$')
plt.plot(x, y[faxial,:], 'bx',     label='solución aproximada por solve_bvp()')
plt.title('Comparación de la solución analítica vs la función solve_bvp() para la carga axial')
plt.xlabel('Eje X (m)')           # titulo del eje X
plt.ylabel('Carga axial (N)')     # titulo del eje Y
plt.legend(loc='upper right')
plt.show()

## bye, bye!!!
