% Calculo de las funciones de forma del elemento placa triangular de Tocher

%% borro memoria, cierro figuras
clear, clc, close all

%% se definen las coordenadas del triangulo a analizar
cx = [ 0 1.0 0.5 ]; cy = [0 0 sqrt(0.75)]; % coordenadas del triangulo

%% se definen las variables simbolicas
syms x y E nu t x1 y1 x2 y2 x3 y3

%% se calculan las matrices PT y A
PT = [ 1 x y x^2 x*y y^2 x^3 (x^2*y+x*y^2) y^3 ];
PT_x = diff(PT,x);
PT_y = diff(PT,y);

A = [ ...
   subs(PT,   {x, y}, {x1, y1})
   subs(PT_x, {x, y}, {x1, y1}) % nodo 1
   subs(PT_y, {x, y}, {x1, y1})

   subs(PT,   {x, y}, {x2, y2})
   subs(PT_x, {x, y}, {x2, y2}) % nodo 2
   subs(PT_y, {x, y}, {x2, y2})

   subs(PT,   {x, y}, {x3, y3})
   subs(PT_x, {x, y}, {x3, y3}) % nodo 3
   subs(PT_y, {x, y}, {x3, y3})
   ]

%% se sustituyen las coordenadas del triangulo en A
A = subs(A,x1,cx(1)); A = subs(A,x2,cx(2)); A = subs(A,x3,cx(3));
A = subs(A,y1,cy(1)); A = subs(A,y2,cy(2)); A = subs(A,y3,cy(3));

%% se calculan las funciones de forma
N = PT/A;  % N = PT*inv(A);

%% se calcula la matriz L
L = [ ...
   diff(PT_x,x)
   diff(PT_y,y)
   diff(PT_x,y) + diff(PT_y,x) ]

%% Reorganizo las funciones de forma en una matriz de 3x3
% donde las filas representan el nodo i
N = [ ...
   N(1)  N(2)  N(3)
   N(4)  N(5)  N(6)
   N(7)  N(8)  N(9) ];

%{
disp('Las funciones de forma son =')
for i = 1:3
   fprintf('N{%d}   = \n',i); pretty(N(i,1));
   fprintf('Nb{%d}  = \n',i); pretty(N(i,2));
   fprintf('Nbb{%d} = \n',i); pretty(N(i,3));
end;
%}

%% Grafico las funciones de forma
LL1 = 0:0.05:1;
LL2 = 0:0.05:1;
[L1,L2] = meshgrid(LL1,LL2);
L3 = 1 - L1 - L2;

L1 = L1(:);
L2 = L2(:);
L3 = round(100*L3(:))/100;
L3(L3<-1e-3) = NaN;

isnanL3 = isnan(L3);
L1(isnanL3) = [];  
L2(isnanL3) = [];  
L3(isnanL3) = [];  

% coordenadas del triangulo
X = L1*cx(1) + L2*cx(2) + L3*cx(3);
Y = L1*cy(1) + L2*cy(2) + L3*cy(3);

TRI = delaunay(X,Y);

for i = 1:3    % contador de nodos
   figure;
   for j = 1:3  % contador de grados de libertad
      subplot(1,3,j);        % Divido el lienzo en 3x1 dibujos
      grid on                % creo la rejilla
      hold on;               % Para que no se sobreescriban los graficos

      % creo una funcion que pueda evaluar numericamente (es como
      % matlabFunction del toolbox simbolico de MATLAB R2010b)
      %  NN = inline(strrep(char(N(i,j)),'*','.*'),'x','y');
      NN = inline(vectorize(char(N(i,j))),'x','y','x1','y1','x2','y2','x3','y3');

      Z = NN(X,Y, cx(1),cy(1), cx(2),cy(2), cx(3),cy(3));
      Z = Z(:);

      trimesh(TRI,X,Y,Z,'LineWidth',2);
      hold on
      trisurf(TRI,X,Y,Z);
      alpha 0.3
      shading interp
      colormap winter
      axis tight
      daspect([1 1 1]);
      view(3);               % vista tridimensional

      xlabel('x', 'FontSize',26); % titulo eje X
      ylabel('y', 'FontSize',26); % titulo eje Y
      switch j % imprimo el titulo
         case 1
            title(sprintf('N_{%d}(x,y)',       i),'FontSize',26);
         case 2
            title(sprintf('(N_b)_{%d}(x,y)',   i),'FontSize',26);
         case 3
            title(sprintf('(N_{bb})_{%d}(x,y)',i),'FontSize',26);
      end;
   end;
end;

return %bye, bye!
