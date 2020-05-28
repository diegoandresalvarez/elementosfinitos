# -*- coding: utf-8 -*-
'''
Programa para calcular la matriz de rigidez K de un elemento como el mostrado
  /|      
  /|________________________________________________
  /|  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^ q(x) -> variable
  /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  /|################################################
  /|                                               o
  /|      E, I, A constante                       /.\
  /|                                            //////

   |-----------------------L-----------------------|
'''

# %%Importación de librerías
import sympy as sp
import numpy as np

# %%Se definen las variables simbólicas y algunas constantes
EI, L, w1, t1, w2, t2 = sp.symbols('EI L w1 t1 w2 t2')
L2 = L**2
L3 = L**3

# %%Se define el vector de desplazamientos nodales a y la matriz K para el
# elemento doblemente empotrado

ae = sp.Matrix([w1, t1, w2, t2])

Keloc = sp.Matrix([[  12*EI/L3,   6*EI/L2,  -12*EI/L3,   6*EI/L2],
                   [   6*EI/L2,   4*EI/L ,   -6*EI/L2,   2*EI/L ],
                   [ -12*EI/L3,  -6*EI/L2,   12*EI/L3,  -6*EI/L2],
                   [   6*EI/L2,   2*EI/L ,   -6*EI/L2,   4*EI/L ]])

# %%Se multiplica Keloc*ae para obtener las 4 ecuaciones correspondientes    
fe = (Keloc*ae)

# La cuarta ecuacion se iguala a cero ya que no hay momentos y se despeja t2
M2 = fe[4-1]     
t2 = sp.solve(M2, t2)

# Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = sp.Matrix([w1, t1, w2, t2[0]])
fe = sp.simplify((Keloc*ae))

# Se elimina el gdl 4 (el de M2) y se recalcula la matriz K:
K = sp.zeros(3)
for i in range(3):
    K[:,i] = fe[:3,:].subs([(w1,int(i==0)),(t1,int(i==1)),(w2,int(i==2))])

# Se imprime la solución
print('La matriz de rigidez para una rótula en el lado derecho es:')
print(f'{sp.pretty(K)}')

###############################################################################

# %%Ahora se hace lo mismo pero con:
#                                  |\
#    ________________________________________________|\
#    ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  |\  q(x) -> variable
#    |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
#    |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
#    ################################################|\
#    o                                               |\
#   /.\                                              |\
#  /////          E, I, A constante                  |\

#    |-----------------------L-----------------------|

# %%Se definen las variables simbólicas y algunas constantes
t1, t2 = sp.symbols('t1 t2')

# %%Se define el vector de desplazamientos nodales a
ae = sp.Matrix([w1, t1, w2, t2])

# %%Se multiplica Keloc*ae para obtener las 4 ecuaciones correspondientes    
fe = (Keloc*ae)

# La segunda ecuacion se iguala a cero ya que no hay momentos y se despeja t1
M1 = fe[2-1]     
t1 = sp.solve(M1, t1)

# Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = sp.Matrix([w1, t1[0], w2, t2])
fe = sp.simplify((Keloc*ae))

# Se elimina el gdl 4 (el de M2) y se recalcula la matriz K:
K = sp.zeros(3)
for i in range(3):
    K[:,i] = fe[[0,2,3],:].subs([(w1,int(i==0)),(w2,int(i==1)),(t2,int(i==2))])

# Se imprime la solución
print('\nLa matriz de rigidez para una rótula en el lado izquierdo es:')
print(f'{sp.pretty(K)}')