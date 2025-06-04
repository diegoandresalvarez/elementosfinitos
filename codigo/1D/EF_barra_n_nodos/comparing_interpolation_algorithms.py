# Interpolation algorithms already implemented in SCIPY
#
# WHO   DATE            WHAT
# DAA   Oct  6, 2020      First algorithm
# DAA   May 30, 2025      Added Akima and Makima
#
# DAA - Diego Andres Alvarez Marin - daalvarez@unal.edu.co

import numpy as np
import matplotlib.pyplot as plt
# %matplotlib qt5
from scipy.interpolate import (
    PchipInterpolator, 
    CubicSpline, 
    Akima1DInterpolator, 
    interp1d, 
    lagrange
)    
# This script shows how to use different interpolation algorithms
# to interpolate a set of points in 1D. 

# The user clicks on the plot to input the data points.
# The algorithms are:
# - PchipInterpolator         : Piecewise Cubic Hermite Interpolating Polynomial
# - CubicSpline               : Cubic spline interpolation
# - Akima1DInterpolator       : Akima's method for 1D interpolation
# - interp1d with 'nearest'   : Nearest neighbor interpolation
# - interp1d with 'linear'    : Linear interpolation
# - interp1d with 'quadratic' : Quadratic interpolation
# - interp1d with kind=5      : Polynomial interpolation of degree 5
# - lagrange                  : Lagrange polynomial interpolation


# Input the data points
n = 10           #number of data points
plt.figure()
plt.title(f'Click with the mouse {n} times. The interpolation will pass through those points')
plt.axis([-5, 5, -5, 5])
pts = np.asarray(plt.ginput(n, timeout=-1, mouse_stop=None))
x, y = pts[:,0], pts[:,1]
plt.close()

algorithms = [
    (PchipInterpolator(  x, y,                       extrapolate=True),         'Pchip'      ),
    (CubicSpline(        x, y, bc_type='not-a-knot', extrapolate=True),         'CubicSpline'),
    (Akima1DInterpolator(x, y, method='akima'),                                 'Akima'      ),
    (Akima1DInterpolator(x, y, method='makima'),                                'Makima'     ),
    (interp1d(           x, y, kind='nearest',       fill_value='extrapolate'), 'nearest'    ),
    (interp1d(           x, y, kind='linear',        fill_value='extrapolate'), 'linear'     ),
    (interp1d(           x, y, kind='quadratic',     fill_value='extrapolate'), 'quadratic'  ),
#   (interp1d(           x, y, kind='cubic',         fill_value='extrapolate'), 'cubic'      ), # same as CubicSpline
    (interp1d(           x, y, kind=5,               fill_value='extrapolate'), '5'          ),
    (lagrange(           x,y),                                                  'Lagrange'   ),
# There are many more functions, as documented in:
# https://docs.scipy.org/doc/scipy/reference/interpolate.html
]
    
plt.figure()
plt.plot(x, y, 'ro')
xi = np.linspace(-10, 10, 501)
for alg in algorithms:
    yi = alg[0](xi)
    if   alg[1] == 'nearest':
        plt.plot(xi, yi, '.',  linewidth=2, label=alg[1])
    elif alg[1] == 'Lagrange':
        plt.plot(xi, yi, 'k-', linewidth=4, label=alg[1])        
    else:
        plt.plot(xi, yi,       linewidth=2, label=alg[1])

plt.legend(loc='best')
plt.ylim([-5, 5])
plt.show()

# More examples here: https://mmas.github.io/interpolation-scipy
