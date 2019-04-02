**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**

[[image:http://www.jeffpalm.com/fox/fox.jpg]]
Fuente: no fuí capaz de encontrarla... esta caricatura es de FOXTROT (http://www.foxtrot.com/). Finalmente observe que falta un \n al final del printf().


=CAPITULO 20. DINAMICA DE ESTRUCTURAS CON ELEMENTOS FINITOS=
==Programa para calcular la matriz de masa consistente de un elemento finito de cercha de dos nodos==
[[code format="matlab"]]
% Defino las variables
syms x x1 x2 rho A L

x2 = x1 + L;

% Defino las funciones de forma
N1 = (x2-x)/L;
N2 = (x-x1)/L;

% Defino la matriz de funciones de forma
N = [ N1 0  N2 0
      0  N1 0  N2 ];

% Se calcula la matriz de masa
M = int(rho*A*N.'*N,x,x1,x2);

disp('M = rho*A*L/6 * ');
disp(M/(rho*A*L/6))
[[code]]

La respuesta del anterior código es:
[[code]]
M = rho*A*L/6 * 
  [ 2, 0, 1, 0]
  [ 0, 2, 0, 1]
  [ 1, 0, 2, 0]
  [ 0, 1, 0, 2]
[[code]]

==Programa para el cálculo de los modos de vibración de cerchas bidimensionales==
Con el siguiente código: [[file:c20_ejemplo_cercha_modos_vibracion.zip]] se calcularon los modos de vibración y se hizo el análisis estático de la cercha mostrada:
[[image:c20_cercha_dinamica.png]]
Se asumió que la densidad del material de la cercha es 7860 kg/m^3 y que todas las barras tienen la misma sección transversal.
[[image:c20_ejemplo_cercha_modos_vibracion.gif]]

==Programa para el cálculo de los modos de vibración de sólidos 3D==
Con el siguiente código: [[file:c20_ejemplo_solido3d_isopar_dinamica.zip]] se calcularon los modos de vibración de un sólido 3D. En la siguiente figura se muestra uno de sus modos de vibración.
[[image:c20_ejemplo_solido3d_mv_01.gif width="800"]]