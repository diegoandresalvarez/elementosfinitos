# -*- coding: utf-8 -*-

# Calculo de las funciones de forma del elemento de viga de Euler-Bernoulli

# %%Importación de librerías

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# %%Definición de variables

x, xi = sp.symbols('x xi')
alpha0, alpha1, alpha2, alpha3 = sp.symbols('alpha0 alpha1 alpha2 alpha3')
w1, w2, dw_dx1, dw_dx2, L, x1, x2 = sp.symbols('w1 w2 dw_dx1 dw_dx2 L x1 x2')

x2 = x1 + L
xm = (x1+x2)/2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% METODO 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('*** Metodo 1 para encontrar las funciones de forma *** \n\n')

# Resulve el sistema de ecuaciones por alpha0, alpha1, alpha2 y alpha3
# aqui cada una de las lineas es igual a cero

r = sp.solve([(alpha0 + alpha1*x1 +   alpha2*x1**2 + alpha3*x1**3 - w1    ),
              (         alpha1    + 2*alpha2*x1 +  3*alpha3*x1**2 - dw_dx1),
              (alpha0 + alpha1*x2 +   alpha2*x2**2 + alpha3*x2**3 - w2    ),
              (         alpha1    + 2*alpha2*x2 +  3*alpha3*x2**2 - dw_dx2)],
             [alpha0, alpha1, alpha2, alpha3])

w = r[alpha0] + r[alpha1]*x + r[alpha2]*x**2 + r[alpha3]*x**3

# Reescribo las funciones de forma en terminos de xi
w = sp.simplify(w.subs([(x, xi*L/2 + xm)]))

# Se factorizan w1, dw_dx1, w2, dw_dx2
w = sp.collect(w, {'w1', 'dw_dx1', 'w2', 'dw_dx2'})

# El programa retorna:
# w = dw_dx1*(L*xi**3/8 - L*xi**2/8 - L*xi/8 + L/8) +
#     dw_dx2*(L*xi**3/8 + L*xi**2/8 - L*xi/8 - L/8) +
#     w1*(xi**3/4 - 3*xi/4 + 1/2) +
#     w2*(-xi**3/4 + 3*xi/4 + 1/2)

# Es decir
N1  = sp.expand(w.subs([(w1, 1), (dw_dx1, 0), (w2, 0), (dw_dx2, 0)])      )
N1b = sp.expand(w.subs([(w1, 0), (dw_dx1, 1), (w2, 0), (dw_dx2, 0)])/(L/2))
N2  = sp.expand(w.subs([(w1, 0), (dw_dx1, 0), (w2, 1), (dw_dx2, 0)])      )
N2b = sp.expand(w.subs([(w1, 0), (dw_dx1, 0), (w2, 0), (dw_dx2, 1)])/(L/2))

# %%Se muestra finalmente w
# Recuerde que
# w(x) = N1(x)w1 + N1b(x)dw_dx1 + N2(x)w2 + N2b(x)dw_dx2
# en la siguiente expresion se pueden ver claramente los terminos
# de N1, N1b, N2 y N2b

print(f'w(xi) = \n{sp.pretty(w)}')
print('es decir:')

print(f'w(xi) = \n{sp.pretty(N1)} \n*w1 + ')
print(f'{sp.pretty(N1b)} \n(L/2)*dw_dx1 + ')
print(f'{sp.pretty(N2)} \n*w2 + ')
print(f'{sp.pretty(N2b)} \n(L/2)*dw_dx2')

print(' ')
print('Siendo las funciones de forma:')
print(f'N_1(xi)  = \n{sp.pretty(N1)}')
print(f'Nb_1(xi) = \n{sp.pretty(N1b)}')
print(f'N_2(xi)  = \n{sp.pretty(N2)}')
print(f'Nb_2(xi) = \n{sp.pretty(N2b)}')

# %%Se grafican las funciones de forma Hermitianas

N1p  = sp.lambdify([xi], N1)
N1bp = sp.lambdify([xi], N1b)
N2p  = sp.lambdify([xi], N2)
N2bp = sp.lambdify([xi], N2b)
xical = np.arange(-1, 1.05, 0.05)  # Puntos de evaluación de las funciones.

fig, ax = plt.subplots()
ax.plot(xical, N1p(xical), linewidth=2, label=r'$N_1(\xi)$')
ax.plot(xical, N1bp(xical), linewidth=2, label=r'$\bar{N}_1(\xi)$')
ax.plot(xical, N2p(xical), linewidth=2, label=r'$N_2(\xi)$')
ax.plot(xical, N2bp(xical), linewidth=2, label=r'$\bar{N}_2(\xi)$')

ax.axis('equal')
ax.set_title('Funciones de forma de la viga de Euler-Bernoulli de dos nodos')
ax.set_xlabel(r'$\xi$')

ax.legend()
ax.grid()
plt.show()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% METODO 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('\n\n*** Metodo 2 para encontrar las funciones de forma *** \n\n')

# %%Definición de variables

xi, L, E, I, q = sp.symbols('xi L E I q')

# Planteo el sistema de ecuaciones y encuentro los coeficientes
x1 = -1
x2 =  1
A = sp.Matrix([[1, x1,   x1**2, x1**3],  # x = [alpha0; alpha1; alpha2; alpha3]
               [0,  1, 2*x1,  3*x1**2],
               [1, x2,   x2**2, x2**3],
               [0,  1, 2*x2,  3*x2**2]])

alpha = A**-1

# %%Defino las funciones de forma

N = sp.zeros(4, 1       )

print('Las funciones de forma son =')
nombre = ['N_1(xi)', 'Nb_1(xi)', 'N_2(xi)', 'Nb_2(xi)']
for i in range(4):
    N[i] = alpha[0, i] + alpha[1, i]*xi + alpha[2, i]*xi**2 + alpha[3, i]*xi**3
    print(f'{nombre[i]}= \n{sp.pretty(N[i])}')

N1  = N[0]
N1b = N[1]
N2  = N[2]
N2b = N[3]

# %%Cálculo la matriz de funciones de forma y su derivada primera y segunda con respecto a xi
NN = sp.expand(sp.Matrix([N1, N1b*L/2, N2, N2b*L/2]).T)
dNN_dxi   = sp.expand(sp.diff(NN, xi))
dNN2_dxi2 = sp.expand(sp.diff(NN, xi, 2))
dNN3_dxi3 = sp.expand(sp.diff(NN, xi, 3))
print(f'\nN = \n{sp.pretty(NN)}')
print(f'\ndNN_dxi = \n{sp.pretty(dNN_dxi)}')
print(f'\ndNN2_dxi2 = \n{sp.pretty(dNN2_dxi2)}')
print(f'\ndNN3_dxi3 = \n{sp.pretty(dNN3_dxi3)}')

# %%Cálculo de la matriz Bb
Bb = sp.simplify(dNN2_dxi2*(4/L**2))
print(f'\nBb = \n{sp.pretty(Bb)}')

# %%Cálculo de la matriz K
K = sp.simplify(sp.integrate(Bb.T*E*I*Bb*L/2, (xi, -1, 1)))
print(f'\nK = (E*I/L^3)* \n{sp.pretty(K/(E*I/L**3))}')

# %%Cálculo del vector de fuerzas nodales equivalentes por fuerzas masicas
fz = q
m  = 0
f  = sp.simplify(sp.integrate(NN.T*fz*L/2 + dNN_dxi.T*m, (xi, -1, 1)))
print(f'\nf = q*L* \n{sp.pretty(f/(q*L))}')

# %%Cálculo de la matriz de masa consistente
rho, A = sp.symbols('rho A')
M = sp.simplify(sp.integrate(rho*A*NN.T*NN*L/2, (xi, -1, 1)))
print(f'\nM = rho*A*L/420* \n{sp.pretty(M/(rho*A*L/420))}')

# # Fin del programa.
