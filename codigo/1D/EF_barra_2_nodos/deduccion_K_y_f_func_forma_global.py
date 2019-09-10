import sympy as sp
from sympy.matrices import Matrix
sp.init_printing(pretty_print=True)

# Defino las variables simbolicas
x, L, E, A, b = sp.symbols('x L E A b')
N1 = sp.Function('N1')
N2 = sp.Function('N2')
N3 = sp.Function('N3')
N4 = sp.Function('N4')

# Defino las funciones de forma
N = Matrix([[N1(x), N2(x), N3(x), N4(x)]]) # matriz de funciones de forma
B = sp.diff(N,x)                           # matriz de deformación
D = E*A                                    # matriz constitutiva

# Matriz de rigidez (ecuación 2.83)
K = sp.simplify(sp.integrate(B.T*D*B, (x, 0, L)))
print('K = '); sp.pprint(K, num_columns=150); print()

# Vector de fuerzas nodales equivalentes (ecuación 2.83)
f = sp.simplify(sp.integrate(N.T*b, (x, 0, L)))
print('f = '); sp.pprint(f)