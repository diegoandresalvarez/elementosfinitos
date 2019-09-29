# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func_EF_T3 import t2ft_T3

# %% CÁLCULO DE UNA VIGA CON ELEMENTOS FINITOS TRIANGULARES PARA TENSION PLANA
# Definición del problema
# Calcule los desplazamientos y las reacciones en los empotramiento, las
# deformaciones y los esfuerzos de la estructura en TENSION PLANA mostrada
# en la figura adjunta

# %% constantes que ayudaran en la lectura del código
X, Y          = 0, 1
NL1, NL2, NL3 = 0, 1, 2

# %% defino las variables/constantes del sólido
Ee   = 200e9        # módulo de elasticidad del sólido (Pa) = 200GPa
nue  = 0.30         # coeficiente de Poisson
te   = 0.10         # espesor del sólido (m)
rhoe = 7850         # densidad (kg/m^3)
g    = 9.81         # aceleración de la gravedad (m/s^2)

# %% Seleccione la malla a emplear
# 1) Malla del ejemplo de la clase
#df   = pd.read_excel('malla_ejemplo.xlsx', sheet_name=None)

# 2) Malla refinada (malla elaborada por David Felipe Cano Perdomo)
df   = pd.read_excel('malla_refinada.xlsx', sheet_name=None)

# %% posición de los nodos:
# xnod: fila=número del nodo, columna=coordenada X=1 o Y=2
xnod = df['xnod'][['x','y']].to_numpy()
nno  = xnod.shape[0]  # número de nodos (número de filas de xnod)

# %% definición de los grados de libertad
ngdl = 2*nno          # número de grados de libertad (dos por nodo)
gdl  = np.reshape(np.arange(ngdl),(nno,2)) # nodos vs grados de libertad

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
LaG = df['LaG'][['NL1','NL2','NL3']].to_numpy() - 1
nef = LaG.shape[0]    # número de EFs (número de filas de LaG)

# %% Relación de cargas puntuales
f   = np.zeros(ngdl)       # vector de fuerzas nodales equivalentes global
cp  = df['carga_punt']
ncp = cp.shape[0]          # número de cargas puntuales
for i in range(ncp):
   f[gdl[cp['nodo'][i]-1, cp['dirección'][i]-1]] = cp['fuerza puntual'][i]

# %% Se dibuja la malla de elementos finitos
plt.figure()
cgx = np.zeros(nef); cgy = np.zeros(nef) # almacena el centro de gravedad
for e in range(nef):
   idx_NL = [NL1, NL2, NL3, NL1]
   plt.plot(xnod[LaG[e, idx_NL], X], xnod[LaG[e, idx_NL], Y], 'b')

   # Calculo la posición del centro de gravedad del triángulo
   cgx[e] = np.mean(xnod[LaG[e,:], X])
   cgy[e] = np.mean(xnod[LaG[e,:], Y])
   plt.text(cgx[e], cgy[e], f'{e+1}', horizontalalignment='center',
                                      verticalalignment='center',   color='b')

plt.plot(xnod[:,X], xnod[:,Y], 'r*')
for i in range(nno):
   plt.text(xnod[i,X], xnod[i,Y], f'{i+1}', color='r')
plt.axis('equal') # tight
plt.title('Malla de elementos finitos')
plt.show()

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K   = np.zeros((ngdl,ngdl)) # matriz de rigidez global como RALA (sparse)                                                         SPARSE
B   = nef * [None]          # contenedor para las matrices de deformación
idx = nef * [None]          # indices asociados a los gdl del EF e

# matriz constitutiva del elemento para TENSION PLANA
De = np.array([[ Ee/(1-nue**2)    , Ee*nue/(1-nue**2),  0              ],
               [ Ee*nue/(1-nue**2), Ee/(1-nue**2)    ,  0              ],
               [ 0                , 0                ,  Ee/(2*(1+nue)) ]])

for e in range(nef):        # ciclo sobre todos los elementos finitos
   # Calculo de la matriz de rigidez del elemento e
   x1 = xnod[LaG[e,NL1],X];              y1 = xnod[LaG[e,NL1],Y]
   x2 = xnod[LaG[e,NL2],X];              y2 = xnod[LaG[e,NL2],Y]
   x3 = xnod[LaG[e,NL3],X];              y3 = xnod[LaG[e,NL3],Y]

   Ae = 0.5*np.linalg.det(np.array([[ 1, x1, y1 ],      #Area del EF e
                                    [ 1, x2, y2 ],
                                    [ 1, x3, y3 ]]))
   if Ae <= 0:
      raise Exception(
         f'La numeración local del EF {e+1} deben especificarse en sentido antihorario.\n')

   # Calculo de la matriz de deformaciones B.
   a1 = x2*y3 - x3*y2;        b1 = y2-y3;        c1 = x3-x2
   a2 = x3*y1 - x1*y3;        b2 = y3-y1;        c2 = x1-x3
   a3 = x1*y2 - x2*y1;        b3 = y1-y2;        c3 = x2-x1

   B[e] = (1/(2*Ae))*np.array([[ b1,  0,   b2,  0,   b3,  0 ],
                               [  0, c1,    0, c2,    0, c3 ],
                               [ c1, b1,   c2, b2,   c3, b3 ]])

   Ke = B[e].T@De@B[e]*te*Ae

   # Calculo del vector de fuerzas nodales equivalentes del elemento e
   # Fuerzas másicas (peso propio)
   fbe = np.array([0, -rhoe*g, 0, -rhoe*g, 0, -rhoe*g]) * Ae*te/3

   fe = fbe # vector de fuerzas nodales equivalentes

   # Ensamblo las contribuciones a las matrices globales
   idx[e] = np.r_[ gdl[LaG[e,NL1],:], gdl[LaG[e,NL2],:], gdl[LaG[e,NL3],:] ]
   K[np.ix_(idx[e],idx[e])] += Ke
   f[np.ix_(idx[e])]        += fe

# %% Muestro la configuración de la matriz K (K es rala)
plt.figure()
plt.spy(K)
plt.title('Los puntos representan los elementos diferentes de cero')
plt.show()

# %% Relación de las cargas superficiales (vector ft)
carga_distr = df['carga_distr']
nlcd = carga_distr.shape[0] # número de lados con carga distribuída
ft   = np.zeros(ngdl)       # fuerzas nodales equivalentes de cargas superficiales                                             SPARSE
for i in range(nlcd):
   e     = carga_distr['elemento'][i] - 1
   lado  = carga_distr['lado'][i]
   carga = carga_distr[['tix','tiy','tjx','tjy']].loc[i].to_numpy()
   fte = t2ft_T3(xnod[LaG[e,:],:], lado, carga, te)

   ft[np.ix_(idx[e])] += fte

# Agrego al vector de fuerzas nodales equivalentes las fuerzas
# superficiales calculadas
f += ft

# %% restricciones y los grados de libertad del desplazamiento conocidos (c)
restric = df['restric']
nres = restric.shape[0]
c    = np.empty(nres, dtype=int)
for i in range(nres):
   c[i] = gdl[restric['nodo'][i]-1, restric['dirección'][i]-1]

# desplazamientos conocidos
ac = restric['desplazamiento'].to_numpy()

# grados de libertad del desplazamiento desconocidos
d = np.setdiff1d(range(ngdl), c)

# %% extraigo las submatrices y especifico las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |  # recuerde que qc=0 (siempre)
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |
Kcc = K[c,:][:,c]; Kcd = K[c,:][:,d]; fd = f[c]
Kdc = K[d,:][:,c]; Kdd = K[d,:][:,d]; fc = f[d]

# %% resuelvo el sistema de ecuaciones
ad = np.linalg.solve(Kdd, fc - Kdc@ac) # calculo desplazamientos desconocidos
qd = Kcc@ac + Kcd@ad - fd              # calculo fuerzas de equilibrio desconocidas

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = np.zeros(ngdl); q = np.zeros(ngdl) # separo la memoria
a[c] = ac;          a[d] = ad          # desplazamientos
q[c] = qd         # q[d] = qc = 0      # fuerzas nodales de equilibrio

# %% imprimo los resultados
tabla_afq = pd.DataFrame(data=np.c_[a, f, q],
                     index=np.arange(ngdl)+1,
                     columns=['a', 'f', 'q'])
print(tabla_afq)

# print('Nodo   Despl_x (m)   Despl_y (m) = ');     [1:nno; reshape(a,2,nno)]'
# print('Nodo Fuerzas nodales equiv. X, Y (N) = '); [1:nno; reshape(f,2,nno)]'
# print('Nodo Fuerzas nodales equil. X, Y (N) = '); [1:nno; reshape(q,2,nno)]'

# %% Dibujo la malla de elementos finitos y las deformada de esta
delta  = np.reshape(a, (nno,2))
escala = 20000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

plt.figure()
for e in range(nef):
   nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]
   plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], 'r')#, lw=0.5)
   plt.plot(xdef[nod_ef, X], xdef[nod_ef, Y], 'b')#, lw=1.0)
plt.gca().set_aspect('equal', adjustable='box')
#legend('Posición original','Posición deformada','Location', 'SouthOutside')
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.title(f'Deformada escalada {escala} veces')
plt.show()
plt.tight_layout()

# %% Se calcula para cada elemento las deformaciones y los esfuerzos
deform = np.zeros((3,nef))
esfuer = np.zeros((3,nef))
for e in range(nef):
   ae = a[idx[e]]               # desplazamientos de los gdl del elemento e
   deform[:,e] = B[e]@ae        # calculo las deformaciones
   esfuer[:,e] = De@deform[:,e] # calculo los esfuerzos

sx = esfuer[0,:];  sy = esfuer[1,:];  txy = esfuer[2,:]
ex = deform[0,:];  ey = deform[1,:];  gxy = deform[2,:]
ez  = -(nue/Ee)*(sx+sy)  # se calculan las deformaciones ez en tensión plana

# %% imprimo y grafico las deformaciones
tabla_exeyezgxy = pd.DataFrame(data=np.c_[ex, ey, ez, gxy],
                     index=np.arange(nef)+1,
                     columns=['ex', 'ey', 'ez', 'gxy'])
print(tabla_exeyezgxy)

fig, ax = plt.subplots()
#subplot(4,1,1); hold on;
for e in range(nef):
   ax.fill(xnod[LaG[e,:],X], xnod[LaG[e,:],Y], ex[e], cmap=plt.cm.jet)
plt.ylabel(r'$\epsilon_x$') #,'FontSize',26)
fig.colorbar(ax=ax)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()

'''
fig, ax = plt.subplots()
im = ax.scatter(x, y, c=c, s=s, cmap=plt.cm.jet)

# Add a colorbar
fig.colorbar(im, ax=ax)

# set the color limits - not necessary here, but good to know how.
im.set_clim(0.0, 1.0)
'''



#https://stackoverflow.com/questions/17660071/python-matplotlib-how-to-use-fill-between-with-a-colormap-to-fill-the-backgroun

'''
subplot(4,1,2); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ey(e))
end;
ylabel('\epsilon_y','FontSize',26); axis equal tight; colorbar;

subplot(4,1,3); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ez(e))
end;
ylabel('\epsilon_z','FontSize',26); axis equal tight; colorbar;

subplot(4,1,4); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),gxy(e))
end;
ylabel('\gamma_{xy}','FontSize',26); axis equal tight; colorbar;



# %% imprimo y grafico los esfuerzos
disp('Esfuerzos (Pa):  (EF,sx,sy,txy) = '); [1:nef; sx; sy; txy]'
figure
subplot(3,1,1); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sx(e))
end;
ylabel('\sigma_x (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(3,1,2); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sy(e))
end;
ylabel('\sigma_y (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(3,1,3); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),txy(e))
end;
ylabel('\tau_{xy} (Pa)','FontSize',26); axis equal tight; colorbar;


# %% Se calculan y grafican para cada elemento los esfuerzos principales y
# %% sus direcciones
# NOTA: esto solo es valido para el caso de TENSION PLANA).
# En caso de DEFORMACION PLANA se deben calcular los valores y vectores
# propios de la matriz de tensiones de Cauchy
#   [dirppales{e}, esfppales{e}] = eig([sx  txy 0    # matriz de esfuerzos
#                                       txy sy  0    # de Cauchy
#                                       0   0   0]);

s1   = (sx+sy)/2 + sqrt(((sx-sy)/2).^2+txy.^2); # esfuerzo normal maximo
s2   = (sx+sy)/2 - sqrt(((sx-sy)/2).^2+txy.^2); # esfuerzo normal minimo
tmax = (s1-s2)/2;                               # esfuerzo cortante maximo
ang  = 0.5*atan2(2*txy, sx-sy); # angulo de inclinacion de s1

# %% Calculo de los esfuerzos de von Mises
s3 = zeros(size(s1));
sv = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2);

# %% imprimo los resultados
disp('Elemento,s1(Pa),s2(Pa),tmax(Pa),angulo(rad) = '); [1:nef; s1; s2; tmax; ang]'
disp('Elemento,Esfuerzos de von Mises (Pa) = '); [1:nef; sv]'

esc = 2; # escala para graficar las flechas
figure
subplot(3,1,1); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s1(e))
end;
# Grafique lineas que indiquen direcciones principales de sigma_1
quiver(cgx,cgy,...     #  En el punto (cgx,cgy) grafique una flecha (linea)
   s1.*cos(ang),s1.*sin(ang),... # indicando la direccion principal de sigma_1
   esc,...                       # con una escala esc
   'k', ...                      # de color negro
  'ShowArrowHead','off',...      # una flecha sin cabeza
  'LineWidth',2,...              # con un ancho de linea 2
  'Marker','.');                 # y en el punto (x,y) poner un punto '.'
quiver(cgx,cgy,...               # La misma flecha ahora en la otra direccion,
   s1.*cos(ang+pi),s1.*sin(ang+pi),...  # es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\sigma_1 (Pa)','FontSize',26); colorbar

subplot(3,1,2); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s2(e))
end;
# Grafique lineas que indiquen direcciones principales de sigma_2
quiver(cgx,cgy,...                         # flecha indicando la direccion
   s2.*cos(ang+pi/2),s2.*sin(ang+pi/2),... # principal de sigma_2
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, s2.*cos(ang-pi/2),s2.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\sigma_2 (Pa)','FontSize',26); colorbar

subplot(3,1,3); hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),tmax(e))
end;
# Grafique lineas que indiquen direcciones principales de tau_max,
quiver(cgx,cgy, tmax.*cos(ang+pi/4),tmax.*sin(ang+pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang-pi/4),tmax.*sin(ang-pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang+3*pi/4),tmax.*sin(ang+3*pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang-3*pi/4),tmax.*sin(ang-3*pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\tau_{max} (Pa)','FontSize',26); colorbar

figure; hold on;
for e in range(nef):
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sv(e))
end;
ylabel('\sigma_v (Pa)','FontSize',26); axis equal tight; colorbar;
title('Esfuerzos de von Mises (Pa)')

# %%
return; # bye, bye!
'''