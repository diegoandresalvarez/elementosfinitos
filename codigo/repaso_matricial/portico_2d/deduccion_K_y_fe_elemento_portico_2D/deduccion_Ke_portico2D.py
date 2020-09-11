# -*- coding: utf-8 -*-

# Programa para deducir la matriz de rigidez de un elemento de pórtico 2D
import sympy as sp

_V_, _M_ = 0, 1

def ensamblar(K, idx, Ke):
    ''' Ensambla la matriz de rigidez local Ke en la matriz de rigidez global K

    Uso:
    ensamblar(K, idx, Ke)

    Parametros de entrada:
    K   -> matriz de rigidez global
    idx -> grados de libertad donde debe proceder el ensamblaje
    Ke  -> matriz de rigidez local
    '''

    # se verifican las dimensiones de las matrices
    nfil, ncol = K.shape
    assert nfil == ncol, "La matriz de rigidez K debe ser cuadrada"

    nfil, ncol = Ke.shape
    assert nfil == ncol == len(idx),\
            "Ke debe ser cuadrada y debe tener el mismo número de filas que idx"

    # se procede con el ensamblaje
    for i in range(nfil):
        for j in range(ncol):
            K[idx[i], idx[j]] += Ke[i,j]

# Se definen las variables simbólicas
w, x, L, EI, EA = sp.symbols('w x L EI EA')
V, M, t, v      = sp.symbols('V M t v', cls=sp.Function)

Kflexion = sp.zeros(4)       # separa memoria para matriz de rigidez global K

for i in range(4):
    sol = sp.dsolve(
           (
               sp.Eq(V(x).diff(x), 0), # se definen las ecuaciones diferenciales
               sp.Eq(M(x).diff(x), V(x)),
               sp.Eq(t(x).diff(x), M(x)/EI),
               sp.Eq(v(x).diff(x), t(x))
           ), 
           [V(x), M(x), t(x), v(x)], 
           ics = {       
               v(0) : int(i==0),    # con sus respectivas condiciones de 
               t(0) : int(i==1),    # frontera  
               v(L) : int(i==2),    
               t(L) : int(i==3)
           }     
    )
  
    Kflexion[:,i] = [[ +sol[_V_].subs(x,0).rhs ],   # Yi  se evaluan las 
                     [ -sol[_M_].subs(x,0).rhs ],   # Mi  reacciones verticales
                     [ -sol[_V_].subs(x,L).rhs ],   # Yj  y los momentos en los
                     [ +sol[_M_].subs(x,L).rhs ]]   # Mj  apoyos    
    
       
Kaxial = EA/L * sp.Matrix([[1, -1],[-1, 1]])

K = sp.zeros(6)
ensamblar(K, [0, 3],       Kaxial)
ensamblar(K, [1, 2, 4, 5], Kflexion) 

sp.pprint(K)