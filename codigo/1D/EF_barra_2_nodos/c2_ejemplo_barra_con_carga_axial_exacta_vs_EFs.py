# -*- coding: utf-8 -*-

# %% definición del problema
# Calcule los desplazamientos y las reacciones en el empotramiento 
# de la barra mostrada
# 
# | b (carga distribuída de magnitud b)
# |->->->->->->->->->->->->->->->->
# |====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
# 1    2    3    4          nno-1  nno
# |<----longitud L de la barra---->|   el area transversal de la barra es A

import numpy as np

# %% defino las variables
nef = 3                       # número de elementos finitos (EF)
nno = nef+1                   # número de nodos
ngdl = nno                    # número de grados de libertad
E   = 200e9     # Pa          # módulo de elasticidad de la barra
A   = (0.01)**2 # m^2         # área transversal de la barra
L   = 2         # m           # longitud de la barra
b   = 1000      # N/m         # fuerza axial aplicada sobre cada EF
P   = 250       # N           # carga nodal al final de la barra

xnod = np.linspace(0, L, nno) # posición de los nodos
Le  = np.diff(xnod)           # longitud de cada EF (= repmat(L/nef, nef, 1))
k   = E*A/Le                  # rigidez de cada EF

# definición de EFs con respecto a nodos
LaG = np.column_stack(([np.arange(1, nno-1 + 1), np.arange(2, nno   + 1)])) - 1

# %% Relación de cargas puntuales
f = np.zeros(ngdl) # vector de fuerzas nodales equivalentes global
f[nno] = P         # relaciono la carga puntual en el nodo "nno"

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K = np.zeros((ngdl,ngdl))   # matriz de rigidez global
for e in range(nef): # ciclo sobre todos los elementos finitos
   idx = LaG[e,:]
   K[np.id_(idx,idx)] += k[e]*np.array([[1. -1.], [-1. 1.]])
   f[idx]             += (b*Le[e]/2)*np.array([1., 1.])

# %% grados de libertad del desplazamiento conocidos y desconocidos
c = 1 - 1
d = np.setdiff1d(np.arange(ngdl), c)

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd | 
#|    | = |         ||    | - |    |     Recuerde que qc=0 (siempre)
#| qc |   | Kdc Kdd || ad |   | fc |

# %% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]; fd = f[c]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = 0               # desplazamientos conocidos (en el gdl 1)

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # calculo desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # calculo fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl)  # separo la memoria
a[c] = ac;          a[d] = ad           # desplazamientos 
q[c] = qd         # q[d] = qc = 0       # fuerzas nodales de equilibrio

# %% calculo las cargas axiales en cada elemento finito
faxial = np.zeros(nef)
for e in range(nef): # ciclo sobre todas los elementos finitos
   Be = np.array([[-1/Le[e], 1/Le[e]]])
   ae = np.array([[a[LaG[e,1]]], [a[LaG[e,2]]]])
   faxial[e] = (E*A)*Be@ae                       # = D*B(e)*a(e)


# %% imprimo los resultados
format short g
disp('Desplazamientos (m) = ');                        a 
disp('Fuerzas nodales equivalentes(N) = ');            f
disp('Fuerzas nodales de equilibrio (N) = ');          q
disp('Cargas axiales en cada elemento finito (N) = '); faxial 

# %% Grafico la solucion analitica y la solucion por el MEF
# 1) grafico los desplazamientos de la barra
u = @(x) (-b*x.^2/2 + (P + b*L)*x)/(E*A); # solucion analitica para el despl.

figure                             # cree un nuevo lienzo
subplot(2,1,1);                    # grafique en la parte superior (1) del lienzo
xx = linspace(0,L,100);            # 100 puntos unif/ distrib. entre 0 y L
plot(xx, u(xx), 'r');              # grafico solucion analitica
hold on;                           # no borre el lienzo 
plot(xnod, a, 'b.-');              # grafico solucion por MEF
title('Comparacion de la solucion analitica con el MEF para el desplazamiento');
xlabel('Eje X (m)')                # titulo del eje X
ylabel('Desplazamiento (m)')       # titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','SouthEast');

# 2) grafico la carga axial de la barra
faxial_exacta = @(x) (P + b*(L-x)); # solucion analitica para la carga axial

subplot(2,1,2);                    # grafique en la parte inferior (2) del lienzo
plot(xx, faxial_exacta(xx), 'r');  # grafico solucion analitica
hold on;                           # no borre el lienzo
for e = 1:nef # ciclo sobre todas los elementos finitos
   plot([xnod(e) xnod(e+1)], [faxial(e) faxial(e)], 'b.-'); # grafico solucion por MEF
end;
title('Comparacion de la solucion analitica con el MEF para la carga axial');
xlabel('Eje X (m)')                # titulo del eje X
ylabel('Carga axial (N)')          # titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','NorthEast');

return; #bye, bye!!!
