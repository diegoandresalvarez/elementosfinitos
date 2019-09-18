# -*- coding: utf-8 -*-

import sympy as sp
from sympy.polys.polyfuncs import interpolate
sp.init_printing(pretty_print=True)

# Nota: el módulo sympy provee unas funciones muy limitadas para la graficación
# Lo mejor será pasar las funciones de forma N_i por sp.lambdify() y graficarlas
# con matplotlib

# %% -------------------------------------------------------------------------
# Funciones de Forma Lagrangianas de dos nodos

# Calculo las funciones de forma
xi = sp.symbols('xi')
N1 = interpolate([(-1,1), (1,0)], xi)
N2 = interpolate([(-1,0), (1,1)], xi)

# Imprimo las funciones de forma
print('\n\nFunciones de Forma Lagrangianas de DOS nodos:')
print(f'\nN1({sp.pretty(xi)}) ='); sp.pprint(N1)
print(f'\nN2({sp.pretty(xi)}) ='); sp.pprint(N2)

# Grafico las funciones de forma
p = sp.plot((N1, (xi, -1, 1)),
            (N2, (xi, -1, 1)), show=False)
p[0].line_color = 'red'
p[0].label = r'$N_1(\xi)$'

p[1].line_color = 'blue'
p[1].label = r'$N_2(\xi)$'

p.xlabel = r'$\xi$'
p.ylabel = ' '
p.title  = 'Funciones de Forma Lagrangianas de DOS nodos'
p.legend = True
p.show()

# %% -------------------------------------------------------------------------
# Funciones de Forma Lagrangianas de tres nodos

# Calculo las funciones de forma
xi = sp.symbols('xi')
N1 = interpolate([(-1,1), (0,0), (1,0)], xi)
N2 = interpolate([(-1,0), (0,1), (1,0)], xi)
N3 = interpolate([(-1,0), (0,0), (1,1)], xi)

# Imprimo las funciones de forma
print('\n\nFunciones de Forma Lagrangianas de TRES nodos:')
print(f'\nN1({sp.pretty(xi)}) ='); sp.pprint(N1)
print(f'\nN2({sp.pretty(xi)}) ='); sp.pprint(N2)
print(f'\nN3({sp.pretty(xi)}) ='); sp.pprint(N3)

# Grafico las funciones de forma
p = sp.plot((N1, (xi, -1, 1)),
            (N2, (xi, -1, 1)),
            (N3, (xi, -1, 1)), show=False)
p[0].line_color = 'red'
p[0].label = r'$N_1(\xi)$'

p[1].line_color = 'blue'
p[1].label = r'$N_2(\xi)$'

p[2].line_color = 'cyan'
p[2].label = r'$N_3(\xi)$'

p.xlabel = r'$\xi$'
p.ylabel = ' '
p.title  = 'Funciones de Forma Lagrangianas de TRES nodos'
p.legend = True
p.show()

# %% -------------------------------------------------------------------------
# Funciones de Forma Lagrangianas de cuatro nodos

# Calculo las funciones de forma
xi = sp.symbols('xi')
N1 = interpolate([(-1,1), (-1/3,0), (1/3,0), (1,0)], xi)
N2 = interpolate([(-1,0), (-1/3,1), (1/3,0), (1,0)], xi)
N3 = interpolate([(-1,0), (-1/3,0), (1/3,1), (1,0)], xi)
N4 = interpolate([(-1,0), (-1/3,0), (1/3,0), (1,1)], xi)

# Imprimo las funciones de forma
print('\n\nFunciones de Forma Lagrangianas de CUATRO nodos:')
print(f'\nN1({sp.pretty(xi)}) ='); sp.pprint(N1)
print(f'\nN2({sp.pretty(xi)}) ='); sp.pprint(N2)
print(f'\nN3({sp.pretty(xi)}) ='); sp.pprint(N3)
print(f'\nN4({sp.pretty(xi)}) ='); sp.pprint(N4)

# Grafico las funciones de forma
p = sp.plot((N1, (xi, -1, 1)),
            (N2, (xi, -1, 1)),
            (N3, (xi, -1, 1)),
            (N4, (xi, -1, 1)), show=False)
p[0].line_color = 'red'
p[0].label = r'$N_1(\xi)$'

p[1].line_color = 'blue'
p[1].label = r'$N_2(\xi)$'

p[2].line_color = 'cyan'
p[2].label = r'$N_3(\xi)$'

p[3].line_color = 'magenta'
p[3].label = r'$N_4(\xi)$'

p.xlabel = r'$\xi$'
p.ylabel = ' '
p.title  = 'Funciones de Forma Lagrangianas de CUATRO nodos'
p.legend = True
p.show()