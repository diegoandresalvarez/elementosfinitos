# -*- coding: utf-8 -*-

# Recuerde escribir primero en la consola %matplotlib qt

# %% El papel del jacobiano en el elemento finito isoparamétrico
import numpy as np
import matplotlib.pyplot as plt
import warnings

# %% Se crea un espacio para hacer clic y definir los nodos del EF
plt.figure()
plt.title('Haciendo clic con el ratón, defina los 16 nodos del EF')
plt.axis([-5, 5, -5, 5])
xnod, ynod = np.array(plt.ginput(16)).T
plt.close()

# %% Espacio normalizado
fig = plt.figure()
gs = fig.add_gridspec(ncols=2, nrows=2)

# se "discretiza" el espacio cada 0.05
xi, eta = np.mgrid[-1:1:41j, -1:1:41j]
n = len(xi)

# se grafica el EF en el espacio normalizado
ax1 = fig.add_subplot(gs[0,0])
for i in range(n):
    h1 = ax1.plot(xi[:, i], eta[:, i], 'b')
    h2 = ax1.plot(xi[i, :], eta[i, :], 'b')
    if i == 0 or i == (n-1):
        plt.setp(h1, linewidth=4)
        plt.setp(h2, linewidth=4)

# se grafican los nodos en el espacio normalizado
xinod =  [-1, -1/3, 1/3,  1,    1,   1, 1, 1/3, -1/3, -1,  -1,   -1, -1/3,  1/3, 1/3, -1/3]
etanod = [-1,   -1,  -1, -1, -1/3, 1/3, 1,   1,    1,  1, 1/3, -1/3, -1/3, -1/3, 1/3,  1/3]
ax1.plot(xinod, etanod, 'ro', markersize=12, linewidth=4)
ax1.axis([-1.1, 1.1, -1.1, 1.1])
ax1.set_aspect('equal', 'box')
ax1.set_xlabel(r'$\xi$')
ax1.set_ylabel(r'$\eta$')

# %% Funciones de forma del elemento lagrangiano plano de 16 nodos (cuadrático)
#
# Numeración local:
#        ^ eta
#        |
#        |
# 10---9---8---7
#  |   |   |   |
# 11--16--15---6
#  |   |   |   |----> xi
# 12--13--14---5
#  |   |   |   |
#  1---2---3---4

Ni1 = (1/16)*(xi - 1)*(1 - 9*xi**2)
Ni2 = (9/16)*(1 - xi**2)*(1 - 3*xi)
Ni3 = (9/16)*(1 - xi**2)*(1 + 3*xi)
Ni4 = (1/16)*(xi + 1)*(9*xi**2 - 1)

Nj1 = (1/16)*(eta - 1)*(1 - 9*eta**2)
Nj2 = (9/16)*(1 - eta**2)*(1 - 3*eta)
Nj3 = (9/16)*(1 - eta**2)*(1 + 3*eta)
Nj4 = (1/16)*(eta + 1)*(9*eta**2 - 1)

N = np.zeros((16, n, n))
N[ 9] = Ni1*Nj4;      N[ 8]=Ni2*Nj4;      N[ 7] = Ni3*Nj4;      N[6] = Ni4*Nj4
N[10] = Ni1*Nj3;      N[15]=Ni2*Nj3;      N[14] = Ni3*Nj3;      N[5] = Ni4*Nj3
N[11] = Ni1*Nj2;      N[12]=Ni2*Nj2;      N[13] = Ni3*Nj2;      N[4] = Ni4*Nj2
N[ 0] = Ni1*Nj1;      N[ 1]=Ni2*Nj1;      N[ 2] = Ni3*Nj1;      N[3] = Ni4*Nj1

# %% Se calcula y grafica el EF en el espacio de geometría real
# se calcula la geometria del EF
x = np.sum(N*xnod[:,np.newaxis,np.newaxis], axis=0) # broadcasting en acción
y = np.sum(N*ynod[:,np.newaxis,np.newaxis], axis=0)

# se grafica el EF en el espacio de geometría real
ax3 = fig.add_subplot(gs[:, 1])
for i in range(n):
    h1 = ax3.plot(x[:, i], y[:, i], 'b')
    h2 = ax3.plot(x[i, :], y[i, :], 'b')
    if i == 0 or i == (n-1):
        plt.setp(h1, linewidth=4)
        plt.setp(h2, linewidth=4)
ax3.plot(xnod, ynod, 'ro', markersize=12, linewidth=4)
ax3.set_aspect('equal', 'box')
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$y$')

# %% derivadas de las funciones de forma
dN_dxi = np.zeros((16, n, n))
dN_dxi[ 0] =      -(-27*xi**2 + 18*xi + 1)*(9*eta**3 - 9*eta**2 - eta + 1)/256
dN_dxi[ 1] =       9*(-9*xi**2 + 2*xi + 3)*(9*eta**3 - 9*eta**2 - eta + 1)/256
dN_dxi[ 2] =        9*(9*xi**2 + 2*xi - 3)*(9*eta**3 - 9*eta**2 - eta + 1)/256
dN_dxi[ 3] =       -(27*xi**2 + 18*xi - 1)*(9*eta**3 - 9*eta**2 - eta + 1)/256
dN_dxi[ 4] =      9*(27*xi**2 + 18*xi - 1)*(3*eta**3 - eta**2 - 3*eta + 1)/256
dN_dxi[ 5] =     -9*(27*xi**2 + 18*xi - 1)*(3*eta**3 + eta**2 - 3*eta - 1)/256
dN_dxi[ 6] =        (27*xi**2 + 18*xi - 1)*(9*eta**3 + 9*eta**2 - eta - 1)/256
dN_dxi[ 7] =       -9*(9*xi**2 + 2*xi - 3)*(9*eta**3 + 9*eta**2 - eta - 1)/256
dN_dxi[ 8] =      -9*(-9*xi**2 + 2*xi + 3)*(9*eta**3 + 9*eta**2 - eta - 1)/256
dN_dxi[ 9] =       (-27*xi**2 + 18*xi + 1)*(9*eta**3 + 9*eta**2 - eta - 1)/256
dN_dxi[10] =    -9*(-27*xi**2 + 18*xi + 1)*(3*eta**3 + eta**2 - 3*eta - 1)/256
dN_dxi[11] =     9*(-27*xi**2 + 18*xi + 1)*(3*eta**3 - eta**2 - 3*eta + 1)/256
dN_dxi[12] =     -81*(-9*xi**2 + 2*xi + 3)*(3*eta**3 - eta**2 - 3*eta + 1)/256
dN_dxi[13] =      -81*(9*xi**2 + 2*xi - 3)*(3*eta**3 - eta**2 - 3*eta + 1)/256
dN_dxi[14] =       81*(9*xi**2 + 2*xi - 3)*(3*eta**3 + eta**2 - 3*eta - 1)/256
dN_dxi[15] =      81*(-9*xi**2 + 2*xi + 3)*(3*eta**3 + eta**2 - 3*eta - 1)/256

dN_deta = np.zeros((16, n, n))
dN_deta[ 0] =      -(-27*eta**2 + 18*eta + 1)*(9*xi**3 - 9*xi**2 - xi + 1)/256
dN_deta[ 1] =     9*(-27*eta**2 + 18*eta + 1)*(3*xi**3 - xi**2 - 3*xi + 1)/256
dN_deta[ 2] =    -9*(-27*eta**2 + 18*eta + 1)*(3*xi**3 + xi**2 - 3*xi - 1)/256
dN_deta[ 3] =       (-27*eta**2 + 18*eta + 1)*(9*xi**3 + 9*xi**2 - xi - 1)/256
dN_deta[ 4] =      -9*(-9*eta**2 + 2*eta + 3)*(9*xi**3 + 9*xi**2 - xi - 1)/256
dN_deta[ 5] =       -9*(9*eta**2 + 2*eta - 3)*(9*xi**3 + 9*xi**2 - xi - 1)/256
dN_deta[ 6] =        (27*eta**2 + 18*eta - 1)*(9*xi**3 + 9*xi**2 - xi - 1)/256
dN_deta[ 7] =     -9*(27*eta**2 + 18*eta - 1)*(3*xi**3 + xi**2 - 3*xi - 1)/256
dN_deta[ 8] =      9*(27*eta**2 + 18*eta - 1)*(3*xi**3 - xi**2 - 3*xi + 1)/256
dN_deta[ 9] =       -(27*eta**2 + 18*eta - 1)*(9*xi**3 - 9*xi**2 - xi + 1)/256
dN_deta[10] =        9*(9*eta**2 + 2*eta - 3)*(9*xi**3 - 9*xi**2 - xi + 1)/256
dN_deta[11] =       9*(-9*eta**2 + 2*eta + 3)*(9*xi**3 - 9*xi**2 - xi + 1)/256
dN_deta[12] =     -81*(-9*eta**2 + 2*eta + 3)*(3*xi**3 - xi**2 - 3*xi + 1)/256
dN_deta[13] =      81*(-9*eta**2 + 2*eta + 3)*(3*xi**3 + xi**2 - 3*xi - 1)/256
dN_deta[14] =       81*(9*eta**2 + 2*eta - 3)*(3*xi**3 + xi**2 - 3*xi - 1)/256
dN_deta[15] =      -81*(9*eta**2 + 2*eta - 3)*(3*xi**3 - xi**2 - 3*xi + 1)/256

# %%Estas derivadas se calcularon con el siguiente código de MATLAB:
'''
import sympy as sp

xi, eta = sp.symbols('xi eta')

Ni1 = sp.nsimplify(sp.interpolate([(-1,1), (-1/3,0), (1/3,0), (1,0)], xi))
Ni2 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,1), (1/3,0), (1,0)], xi))
Ni3 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,0), (1/3,1), (1,0)], xi))
Ni4 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,0), (1/3,0), (1,1)], xi))

Nj1 = sp.nsimplify(sp.interpolate([(-1,1), (-1/3,0), (1/3,0), (1,0)], eta))
Nj2 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,1), (1/3,0), (1,0)], eta))
Nj3 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,0), (1/3,1), (1,0)], eta))
Nj4 = sp.nsimplify(sp.interpolate([(-1,0), (-1/3,0), (1/3,0), (1,1)], eta))

N = 16 * [None]

N[ 9] = Ni1 * Nj4; N[ 8] = Ni2 * Nj4; N[ 7] = Ni3 * Nj4; N[6] = Ni4 * Nj4; 
N[10] = Ni1 * Nj3; N[15] = Ni2 * Nj3; N[14] = Ni3 * Nj3; N[5] = Ni4 * Nj3;
N[11] = Ni1 * Nj2; N[12] = Ni2 * Nj2; N[13] = Ni3 * Nj2; N[4] = Ni4 * Nj2;
N[ 0] = Ni1 * Nj1; N[ 1] = Ni2 * Nj1; N[ 2] = Ni3 * Nj1; N[3] = Ni4 * Nj1;

dN_dxi  = 16 * [None]
dN_deta = 16 * [None]

for i in range(16):
    dN_dxi[i] = sp.simplify(sp.diff(N[i], xi))
    print(f'dN_dxi[{i:2}] = {str(dN_dxi[i]):>65}')

for i in range(16):
    dN_deta[i] = sp.simplify(sp.diff(N[i], eta))
    print(f'dN_deta[{i:2}] = {str(dN_deta[i]):>64}')
'''
#%% Calculo el determinante del Jacobiano
dx_dxi  = np.sum(dN_dxi *xnod[:,np.newaxis,np.newaxis], axis=0) # broadcasting en acción
dy_dxi  = np.sum(dN_dxi *ynod[:,np.newaxis,np.newaxis], axis=0) 
dx_deta = np.sum(dN_deta*xnod[:,np.newaxis,np.newaxis], axis=0) 
dy_deta = np.sum(dN_deta*ynod[:,np.newaxis,np.newaxis], axis=0) 

# J = [ dx_dxi   dy_dxi
#       dx_deta  dy_deta ]

detJ = dx_dxi*dy_deta - dx_deta*dy_dxi

# Se calcula el Jacobian ratio
if np.max(detJ) > 0:
    JR = np.max(detJ) / np.min(detJ)
else:
    JR = np.min(detJ) / np.max(detJ)

if JR > 0:
    print(f'Jacobian ratio (JR) = {JR}')
print(f'min(detJ) = {np.min(detJ)}')
print(f'max(detJ) = {np.max(detJ)}')

if np.max(detJ) < 0:
    print('El EF está numerado en sentido contrario al especificado')

if JR < 0 or JR > 40:
    print('Esta forma no es adecuada para un EF.')

# %% se grafica el jacobiano
ax2 = fig.add_subplot(gs[1,0])
h = ax2.pcolor(xi, eta, detJ, cmap='jet',shading='auto')
tickscb = np.linspace(detJ.min(), detJ.max(), 10)
fig.colorbar(h, ax=ax2, ticks=tickscb), 
ax2.set_title('Determinante de J')
ax2.set_xlabel(r'$\xi$')
ax2.set_ylabel(r'$\eta$')

if JR < 0 or JR > 40:
    ax3.set_title(f'Esta forma no es adecuada para un EF.\nJacobian ratio = {JR:.3}')
else:
    ax3.set_title(f'Jacobian ratio = {JR:.3}')     
ax2.set_aspect('equal', 'box')
warnings.filterwarnings("ignore")
ax2.contour(xi, eta, detJ, levels= [0] , linewidths=4)
warnings.filterwarnings("default")
plt.show()
