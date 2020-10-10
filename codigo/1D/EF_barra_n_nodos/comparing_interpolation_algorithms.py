# Interpolation algorithms already implemented in SCIPY
#
# WHO   DATE            WHAT
# DAA   Oct 6, 2020      First algorithm
#
# DAA - Diego Andres Alvarez Marin - diegoandresalvarez@gmx.net

import numpy as np
import matplotlib.pyplot as plt
# %matplotlib qt5
import scipy.interpolate as interpolate
from numpy.polynomial.polynomial import Polynomial


# Input the data points
n = 10           #number of data points

plt.figure()
plt.title(f'Click with the mouse {n} times. The interpolation will pass through those points')
plt.axis([-5, 5, -5, 5])
pts = np.asarray(plt.ginput(n, timeout=-1, mouse_stop=None))
x, y = pts[:,0], pts[:,1]
plt.close()

algorithms = [
(interpolate.interp1d(x, y, kind='nearest',   fill_value='extrapolate'), 'nearest'  ),    
(interpolate.interp1d(x, y, kind='linear',    fill_value='extrapolate'), 'linear'   ),
(interpolate.interp1d(x, y, kind='quadratic', fill_value='extrapolate'), 'quadratic'),
(interpolate.interp1d(x, y, kind='cubic',     fill_value='extrapolate'), 'cubic'    ),
(interpolate.interp1d(x, y, kind=5,           fill_value='extrapolate'), '5'        ),
(interpolate.lagrange(x,y),                                              'Lagrange' ),
# There are many more functions, as documented in :
# https://docs.scipy.org/doc/scipy/reference/interpolate.html
]
    
plt.figure()
plt.plot(x, y, 'ro')
xi = np.linspace(-10, 10, 501)
for alg in algorithms:
    yi = alg[0](xi)
    if alg[1] == 'nearest':
        plt.plot(xi, yi, '.', linewidth=2, label=alg[1])
    if alg[1] == 'Lagrange':
        plt.plot(xi, yi, 'k-', linewidth=4, label=alg[1])        
    else:
        plt.plot(xi, yi,      linewidth=2, label=alg[1])

plt.legend(loc='best')
plt.ylim([-5, 5])
plt.show()

# más ejemplos aquí: https://mmas.github.io/interpolation-scipy

