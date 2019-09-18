# Integration using a Gauss-Legendre quadrature:

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import legendre
import matplotlib.pyplot as plt

m = 4

xi, w = leggauss(m) # esta función reemplaza un posible "gausslegendre_quad.py"

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

a, b = 0, 0.8
f = lambda x : 0.2 + 25*x - 200*x**2 +675*x**3 - 900*x**4 + 400*x**5
print('Error =', abs(((b-a)/2)*np.sum(w*f((b+a)/2 + (b-a)*xi/2)) - 3076/1875))

a, b = 0, np.pi/2
f = lambda x : np.sin(x)
print('Error =', abs(((b-a)/2)*np.sum(w*f((b+a)/2 + (b-a)*xi/2)) - 1))

# NOTA: se recomienda utilizar la función scipy.integrate.fixed_quad
# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.fixed_quad.html#scipy.integrate.fixed_quad
