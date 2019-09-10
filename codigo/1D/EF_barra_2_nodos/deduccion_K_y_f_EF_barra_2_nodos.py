import sympy as sp
from sympy.matrices import Matrix

# Defino las variables simbolicas
x, x1, x2, E, A, L, b = sp.symbols('x x1 x2 E A L b')

# Defino las funciones de forma
x2 = x1 + L
N1 = (x2-x)/L;                 N2 = (x-x1)/L

N = Matrix([[N1, N2]])                        # matriz de funciones de forma
B = Matrix([[sp.diff(N1,x), sp.diff(N2,x)]])  # matriz de deformación
D = E*A                                       # matriz constitutiva

# Matriz de rigidez (ecuación 2.83)
K = sp.simplify(sp.integrate(B.T*D*B, (x, x1, x2)))
print('K = '); sp.pprint(K); print()

# Vector de fuerzas nodales equivalentes (ecuación 2.83)
f = sp.simplify(sp.integrate(N.T*b, (x, x1, x2)))
print('f = '); sp.pprint(f)