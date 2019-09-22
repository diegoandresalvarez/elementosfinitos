function ft = t2ft_T3(xnod, lado, carga, espesor)
% Esta función convierte las fuerzas superficiales aplicadas a un elemento
% finito triangular de 3 nodos a sus correspondientes cargas nodales
% equivalentes ft
%
% xnod = [ x1e y1e
%          x2e y2e
%          x3e y3e ];
%
% lado = 12, 23, 31
%
% carga = [ t1x t1y t2x t2y ];   % si la carga se aplica sobre el lado 12
%         [ t2x t2y t3x t3y ];   % si la carga se aplica sobre el lado 23
%         [ t3x t3y t1x t1y ];   % si la carga se aplica sobre el lado 31

%% Se definen algunas constantes
X = 1; Y = 2;

switch lado
   case 12
      %% Fuerzas sobre el lado 12
      % Se calcula la longitud del lado 12
      L12 = sqrt((xnod(1,X) - xnod(2,X))^2 + (xnod(1,Y) - xnod(2,Y))^2);
      
      % Fuerzas distribuidas aplicadas en los nodos 1 y 2 locales
      t1x = carga(1); t1y = carga(2); t2x = carga(3); t2y = carga(4);
      ft = espesor*[ ...
         (L12*(2*t1x + t2x))/6
         (L12*(2*t1y + t2y))/6
         (L12*(t1x + 2*t2x))/6
         (L12*(t1y + 2*t2y))/6
         0
         0 ];
      
   case 23
      %% Fuerzas sobre el lado 23
      % Se calcula la longitud del lado 23
      L23 = sqrt((xnod(2,X) - xnod(3,X))^2 + (xnod(2,Y) - xnod(3,Y))^2);
      
      % Fuerzas distribuidas aplicadas en los nodos 2 y 3 locales
      t2x = carga(1); t2y = carga(2); t3x = carga(3); t3y = carga(4);
      ft = espesor*[ ...
         0
         0
         (L23*(2*t2x + t3x))/6
         (L23*(2*t2y + t3y))/6
         (L23*(t2x + 2*t3x))/6
         (L23*(t2y + 2*t3y))/6 ];
      
   case 31
      %% Fuerzas sobre el lado 31
      % Se calcula la longitud del lado 31
      L31 = sqrt((xnod(3,X) - xnod(1,X))^2 + (xnod(3,Y) - xnod(1,Y))^2);
      
      % Fuerzas distribuidas aplicadas en los nodos 3 y 1 locales
      t3x = carga(1); t3y = carga(2); t1x = carga(3); t1y = carga(4);
      ft = espesor*[ ...
         (L31*(2*t1x + t3x))/6
         (L31*(2*t1y + t3y))/6
         0
         0
         (L31*(t1x + 2*t3x))/6
         (L31*(t1y + 2*t3y))/6 ];
   otherwise
      error('Unicamente se permiten los lados 12, 23 o 31');
end

%%
return;


%% NOTA: Los vectores se calcularon con el siguiente programa:
%{
%% Elemento triangular de 3 nodos
clear
clc
syms N1 N2 N3 
syms t1x t1y t2x t2y t3x t3y 
syms xi L12 L23 L31

N = [ N1 0   N2 0   N3 0  
      0  N1  0  N2  0  N3 ];
 
t = [ t1x t1y t2x t2y t3x t3y ].';

res = N.'*N*t;
 
res12 = res;
res12 = subs(res12, N1, (1-xi)/2);
res12 = subs(res12, N2, (1+xi)/2);
res12 = subs(res12, N3, 0);
dx_dxi = L12/2;
f12 = int(res12*dx_dxi, xi, -1 ,1)

res23 = res;
res23 = subs(res23, N2, (1-xi)/2);
res23 = subs(res23, N3, (1+xi)/2);
res23 = subs(res23, N1, 0);
dx_dxi = L23/2;
f23 = int(res23*dx_dxi, xi, -1 ,1)

res31 = res;
res31 = subs(res31, N3, (1-xi)/2);
res31 = subs(res31, N1, (1+xi)/2);
res31 = subs(res31, N2, 0);
dx_dxi = L31/2;
f31 = int(res31*dx_dxi, xi, -1 ,1)
%}