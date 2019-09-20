# -*- coding: utf-8 -*-

# La función
from scipy.integrate import fixed_quad
# realiza integraciones de funciones utilizando cuadraturas de Gauss-Legendre,
# de forma similar a como lo hace el programa de MATLAB "cuadratura_poly_sin_GL.m".
# La idea de este programa es mostrar como se utiliza dicha función.

import numpy as np
import matplotlib.pyplot as plt


# Integración de:
#
f   = lambda x : 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
#
# entre 0 y 0.8 usando cuadraturas de Gauss-Legendre

a   = 0                    # límites de integración
b   = 0.8        
sol = 3076/1875            # solución exacta
err = np.zeros(10)         # separo la memoria
for m in range(10):        # varío el número de puntos de la cuadratura
   err[m] = abs(fixed_quad(f, a, b, n=m+1)[0] - sol)

plt.figure()               # creo un lienzo
plt.plot(range(1,11), err) # grafico el error
plt.xlabel('Número de puntos en la cuadratura')
plt.ylabel('Error')
plt.title('Cuadratura de Gauss Legendre')
plt.grid()                 # pongo la rejilla
plt.show()

# Integración de:
#
f   = lambda x : np.sin(x)
#
# entre 0 y pi/2 usando cuadraturas de Gauss-Legendre

a   = 0                    # límites de integración
b   = np.pi/2
sol = 1                    # solución exacta
err = np.zeros(10)         # separo la memoria
for m in range(10):        # varío el número de puntos de la cuadratura
   err[m] = abs(fixed_quad(f, a, b, n=m+1)[0] - sol)

plt.figure()               # creo un lienzo
plt.semilogy(range(1,11), err) # escala logarítmica para apreciar mejor el error
plt.xlabel('Número de puntos en la cuadratura')
plt.ylabel('Error')
plt.title('Cuadratura de Gauss Legendre')
plt.grid()                 # pongo la rejilla
plt.show()
