function [Te,Ke] = Te_Keloc_RR(A, E, x1, y1, x2, y2)
% Calcula de la matriz de transformacion T y de la matriz de rigidez 
% local Ke de una barra con 
% rotula a la izquierda y con 
% rotula a la derecha

%% Se calcula la longitud de la barra
L = sqrt((x2-x1)^2 + (y2-y1)^2);  

%% Se calcula la matriz T
c = (x2-x1)/L;   s = (y2-y1)/L;  % sin y cos de la inclinacion

Te = [ c  s  0  0
      -s  c  0  0
       0  0  c  s
       0  0 -s  c ];

%% Se calcula la matriz K local
k = A*E/L;

Ke = [ k   0  -k   0       
       0   0   0   0
      -k   0   k   0
       0   0   0   0 ];      

return