# -*- coding: utf-8 -*-

# Integración utilizando cuadraturas de Gauss-Legendre
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import legendre
import matplotlib.pyplot as plt

# %% orden de la cuadratura a utilizar
m = 4

# %% Cálculo de los polinomios de Legendre usando la recursión de Bonnet
P = [None]*(m+1)
P[0] = np.array([1])
P[1] = np.array([1, 0])
for n in range(2, m+1):
   P[n] = ((2*n-1)*np.r_[P[n-1], 0] - (n-1)*np.r_[0, 0, P[n-2]])/n

# %% Raíces del polinomio de Legendre
xi = np.sort(np.roots(P[m]))

# %% Cálculo de los pesos
dPm_dxi = np.polyder(P[m])
w = 2/((1 - xi**2)*np.polyval(dPm_dxi,xi)**2)

# %% Alternativamente, se puede utilizar la función leggauss()
xi2, w2 = leggauss(m)

# %% Se grafican los polinomios de Legendre
x = np.linspace(-1, 1, 100)
for n in range(m+1):
    y = legendre(n)(x)
    plt.plot(x, y, label=f'n={n}')

plt.xlim(-1.0, 1.0)
plt.ylim(-1.1, 1.1)
plt.legend()
plt.grid()
plt.title(f'Primeros {m+1} polinomios de Legendre')
plt.plot(xi, np.zeros_like(xi), 'ro')
plt.show()

# %% Se integran dos ejemplos:
a, b = 0, 0.8
f    = lambda x : 0.2 + 25*x - 200*x**2 +675*x**3 - 900*x**4 + 400*x**5
sol  = 3076/1875
print('Error =', abs(((b-a)/2)*np.sum(w*f((b+a)/2 + (b-a)*xi/2)) - sol))

a, b = 0, np.pi/2
f    = lambda x : np.sin(x)
sol  = 1
print('Error =', abs(((b-a)/2)*np.sum(w*f((b+a)/2 + (b-a)*xi/2)) - sol))

# NOTA: se recomienda utilizar la función scipy.integrate.fixed_quad()
# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.fixed_quad.html#scipy.integrate.fixed_quad
