# -*- coding: utf-8 -*-
# Imposición de las deformaciones angulares gxz
# Se realiza para el elemento finito de Timoshenko de tres nodos

# %%Importación  de librerías
from sympy import Matrix, symbols, interpolate, sqrt, integrate, pretty, nsimplify

# %%Definción de variables
G, Aast, L, xi, w1, w2, w3, t1, t2, t3 = symbols('G Aast L xi w1 w2 w3 t1 t2 t3')

# La matriz Bs "normal"
Bs = lambda xi: (2/L)*Matrix([ xi - 1/2, -(L*xi*(xi - 1))/4, -2*xi,
                              (L*(xi**2 - 1))/2, xi + 1/2, -(L*xi*(xi + 1))/4]).T

# Las funciones de interpolacion de las deformaciones gxz                        
Ng1 = interpolate([(-1/sqrt(3), 1), (1/sqrt(3), 0)], xi) 
Ng2 = interpolate([(-1/sqrt(3), 0), (1/sqrt(3), 1)], xi)

# Se calcula la matriz Bc_sustitutiva
Bs_sustitutiva = nsimplify(Ng1*Bs(-1/sqrt(3)) + Ng2*Bs(+1/sqrt(3)))

# y la correspondiente matriz de rigidez
Ks = integrate(Bs_sustitutiva.T*G*Aast*Bs_sustitutiva*L/2, (xi,-1,1))

print(f'Ks = ((G*Aast)/(9*L)) * \n{pretty(Ks/(G*Aast/(9*L)))}')

# Este Ks coincide con aquel integrado con dos puntos de GL.

# Fin del programa