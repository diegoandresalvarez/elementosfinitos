clear,clc        % borra la memoria y la pantalla

%% Programa para calcular los modos de energía nula del sólido

a = 1;        % ancho elemento (m)
b = 1;        % altura elemento (m)
E  = 200;     % modulo de elasticidad del elemento (GPa)
nu = 0.33;    % coeficiente de Poisson
t  = 0.10;    % espesor del elemento (m)

syms r s

%% Funciones de forma del rectangulo
N1 = (1-r/a)*(1-s/b)/4;
N2 = (1+r/a)*(1-s/b)/4;
N3 = (1+r/a)*(1+s/b)/4;
N4 = (1-r/a)*(1+s/b)/4;

%% matriz constitutiva del elemento para TENSION PLANA
D = [ E/(1-nu^2)     E*nu/(1-nu^2)  0
      E*nu/(1-nu^2)  E/(1-nu^2)     0
      0              0              E/(2*(1+nu)) ];

  %% matriz de deformaciones
B1 = [diff(N1,r)   0         
      0            diff(N1,s)
      diff(N1,s)   diff(N1,r)];
      
B2 = [diff(N2,r)   0
      0            diff(N2,s)
      diff(N2,s)   diff(N2,r)];
      
B3 = [diff(N3,r)   0
      0            diff(N3,s)
      diff(N3,s)   diff(N3,r)];
      
B4 = [diff(N4,r)   0
      0            diff(N4,s)
      diff(N4,s)   diff(N4,r)];

% se ensambla la matriz de deformaciones
B = [B1 B2 B3 B4];

%% realizo la integracion
K = double(int(int(B.'*D*B*t, r, -a,a), s,-b,b));

%% Se calculan los valores y vectores propios de la matriz K
[evec,eval] = eig(K);

%% Coordenadas del elemento
xnod = [ ...
   -1 -1
    1 -1
    1  1
   -1  1 ];

%% Se imprimen los vectores propios (recuerde que los modos de energía nula
%% son aquellos para los cuales los valores propios son cero
figure
modo = cell(8,1);
for i = 1:8
   modo{i} = reshape(evec(:,i),2,4)' + xnod;
   subplot(2,4,i)
   plot(xnod([1 2 3 4 1],1),xnod([1 2 3 4 1],2),'b');
   hold on
   plot(modo{i}([1 2 3 4 1],1),modo{i}([1 2 3 4 1],2),'r');
   axis equal
   axis([-2 2 -2 2]);
   title(sprintf('\\lambda_%d = %d',i,eval(i,i)),'FontSize',16);
end
