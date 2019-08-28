# -*- coding: utf-8 -*-

import sympy as sp

# Se definen las variables simbólicas
E, A, L, u3, u4, P = sp.symbols('E A L u3 u4 P')

# Se calculan las rigideces para cada barra
k1 = E*A/L; k2 = E*A/L; k3 = 2*E*A/L

# Se define la matriz de rigidez
K = sp.Matrix([
    [ k1,   0 ,  -k1      ,   0  ],
    [ 0 ,   k2,  -k2      ,   0  ],
    [-k1,  -k2,   k1+k2+k3,  -k3 ],
    [ 0 ,   0 ,  -k3      ,   k3 ]])

# Se define el vector de desplazamientos nodales, teniendo en cuenta que
# u1=0 y u2=0
a = sp.Matrix([0, 0, u3, u4])  # SymPy lo toma como un vector columna

# Se define el vector de fuerzas nodales de equilibrio
f = sp.Matrix([ 0, 0, 0, P])

# Se definen g.d.l. conocidos y desconocidos asociados a los desplazamientos
c = sp.Matrix([1, 2]) - sp.ones(2,1);       d = sp.Matrix([3, 4]) - sp.ones(2,1)
#c = sp.Matrix([1, 2]).applyfunc(lambda x : x-1)

# Se descomponen los vectores a, f y la matriz K
Kcc = K.extract(c,c);      Kcd = K.extract(c,d)
Kdc = K.extract(d,c);      Kdd = K.extract(d,d)
ac  = a.extract(c,[0]);  # ad  = a.extract(d,[0])
fc  = f.extract(d,[0]);    fd  = f.extract(c,[0])

# Se calculan los vectores ad y fd
ad = Kdd.solve(fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad - fd

sp.pprint(ad)
sp.pprint(qd)