# -*- coding: utf-8 -*-

# Cálculo de las funciones de forma del elemento de viga de Timoshenko lineal.

# %%Importación de librerías
from sympy import symbols, Matrix, diff, integrate, pretty, sqrt, nsimplify, simplify

# %%Defino las variables.
x, xi, E, I, G, Aast, L = symbols('x xi E I G Aast L')
w1, t1, w2, t2, fz, m = symbols('w1 t1 w2 t2 fz m')

# %%Defino las funciones de forma.
N1 = (1-xi)/2
N2 = (1+xi)/2

# %%Defino las matrices de deformación.
Bb = Matrix([ 0,                 (2/L)*diff(N1,xi), 0,                 (2/L)*diff(N2,xi) ]).T
Bs = Matrix([ (2/L)*diff(N1,xi), -N1,               (2/L)*diff(N2,xi), -N2               ]).T

print(f'\nBb = \n{pretty(Bb)}')
print(f'\nBb = \n{pretty(Bs)}')
print('\n\n\n')

# %%Defino las matrices de rigidez
Kb = integrate(Bb.T*E*I*Bb*L/2, (xi, -1, 1))
Ks = integrate(Bs.T*G*Aast*Bs*L/2, (xi, -1, 1))

print('\nSolución exacta: ')
print(f'\nKb = (E*I/L) * \n{pretty(Kb/(E*I/L))}')
print(f'\nKs = (G*Aast/L) * \n{pretty(Ks/(G*Aast/L))}')
print('\n\n\n\n\n')

# %%Evalúo las matrices con una cuadratura de Gauss-Legendre de orden 1.
w1 = 2; xi1 = 0
Kb = (Bb.T*E*I*Bb*L/2).subs(xi, xi1)*w1
Ks = (Bs.T*G*Aast*Bs*L/2).subs(xi, xi1)*w1

print('Integrando con GL de orden 1 : ')
print(f'\nKb = (E*I/L) * \n{pretty(Kb/(E*I/L))}')
print(f'\nKs = (G*Aast/L) * \n{pretty(Ks/(G*Aast/L))}')
print('\n\n\n\n\n')

# %%Evalúo las matrices con una cuadratura de Gauss-Legendre de orden 2.
w1 = 1; xi1 = -sqrt(1/3)
w2 = 1; xi2 = +sqrt(1/3)
Kb = nsimplify((Bb.T*E*I*Bb*L/2).subs(xi, xi1)*w1 +
               (Bb.T*E*I*Bb*L/2).subs(xi, xi2)*w2)

Ks = nsimplify((Bs.T*G*Aast*Bs*L/2).subs(xi, xi1)*w1 +
               (Bs.T*G*Aast*Bs*L/2).subs(xi, xi2)*w2)

print('Integrando con GL de orden 2 : ')
print(f'\nKb = (E*I/L) * \n{pretty(Kb/(E*I/L))}')
print(f'\nKs = (G*Aast/L) * \n{pretty(Ks/(G*Aast/L))}')
print('\n\n')

# %%Evalúo la curvatura.
kappa = simplify((2/L)*(diff(N1, xi)*t1 + diff(N2, xi)*t2))
print(f'kappa = \n{pretty(kappa)}')
print('\n\n')

# %%Evalúo gamma_xz.
gxz = simplify((2/L)*(diff(N1, xi)*w1 + diff(N2, xi)*w2) - (N1*t1 + N2*t2))
print(f'gxz =  \n{pretty(gxz)}')
print('\n\n')

# %%El vector de fuerzas nodales equivalentes fe asociados a una carga
#   distribuída q constante y un momento distribuído m tambien uniforme.
f1 = integrate(N1*Matrix([fz, m]) * L/2, (xi, -1, 1))
f2 = integrate(N2*Matrix([fz, m]) * L/2, (xi, -1, 1))
fe = Matrix([f1, f2])
print('Para una carga distribuida de magnitud q constante y momento'
      + ' distribuido m uniforme')
print(f'\nfe = \n{pretty(fe)}')
print('\n\n')

# %%Calculo el vector de fuerzas nodales equivalentes correspondientes a una
#   carga distribuida trapezoidal y sin momentos distribuidos
q1, q2 = symbols('q1 q2')
fz = N1*q1 + N2*q2
m = 0
f1 = integrate(N1*Matrix([fz, m])*L/2, (xi, -1, 1))
f2 = integrate(N2*Matrix([fz, m])*L/2, (xi, -1, 1))
fe = Matrix([f1, f2])
print('Para una carga distribuida trapezoidal y sin momentos distribuidos')
print(f'\nfe = \n{pretty(fe)}')

# Fin del programa
