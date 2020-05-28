% Calculo de las funciones de forma del elemento placa triangular de Tocher

%% borro memoria, cierro figuras
clear, clc, close all

%% se definen las coordenadas del triangulo a analizar
cx = [ 0 1.0 0.5 ]; cy = [0 0 sqrt(0.75)]; % coordenadas del triangulo

%% se definen las variables simbolicas
syms x y E nu t x1 y1 x2 y2 x3 y3

%% se calculan las matrices PT y A
pT = [ 1 x y x^2 x*y y^2 x^3 (x^2*y+x*y^2) y^3 ];
pT_x = diff(pT,x);
pT_y = diff(pT,y);

A = [ ...
   subs(pT,   {x, y}, {x1, y1})
   subs(pT_x, {x, y}, {x1, y1}) % nodo 1
   subs(pT_y, {x, y}, {x1, y1})

   subs(pT,   {x, y}, {x2, y2})
   subs(pT_x, {x, y}, {x2, y2}) % nodo 2
   subs(pT_y, {x, y}, {x2, y2})

   subs(pT,   {x, y}, {x3, y3})
   subs(pT_x, {x, y}, {x3, y3}) % nodo 3
   subs(pT_y, {x, y}, {x3, y3})
   ]

%% se sustituyen las coordenadas del triangulo en A
A = subs(A, {x1, x2, x3}, cx);
A = subs(A, {y1, y2, y3}, cy);

%% se calculan las funciones de forma
N = pT/A;  % N = pT*inv(A);
Nvec = N;

%% se calcula la matriz L
L = [ ...
     diff(pT,x,2)
     diff(pT,y,2)
   2*diff(pT_x,y) ]

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
      subplot(1,3,j);        % divido el lienzo en 3x1 dibujos
      grid on                % creo la rejilla
      hold on                % para que no se sobreescriban los graficos

      NN = matlabFunction(N(i,j), 'Vars', {'x','y','x1','y1','x2','y2','x3','y3'});      

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

      xlabel('$x$', 'FontSize',26, 'interpreter','latex'); % titulo eje X
      ylabel('$y$', 'FontSize',26, 'interpreter','latex'); % titulo eje Y
      switch j % imprimo el titulo
         case 1
            title(sprintf('$N_{%d}(x,y)$',       i), ...
                'FontSize',26, 'interpreter','latex');
         case 2
            title(['$\overline{N}_' num2str(i) '(x,y)$'], ...
                'FontSize',26, 'interpreter','latex');
         case 3
            title(['$\overline{\overline N}_' num2str(i) '(x,y)$'], ...
                'FontSize',26, 'interpreter','latex');
      end
   end
end

%% Se calcula una matriz QQ que permite calcular las fuerzas cortantes
d3N_dx3   = diff(Nvec,x,3); 
d3N_dx2dy = diff(Nvec,x,x,y);
d3N_dxdy2 = diff(Nvec,x,y,y);
d3N_dy3   = diff(Nvec,y,3);

% Se calcula la matriz QQ, de modo que Q = -D*QQ*a
QQ = simplify([ d3N_dx3 + d3N_dxdy2
                d3N_dy3 + d3N_dx2dy ])

return %bye, bye!
