% Calculo de las funciones de forma del elemento placa MZC

clear, clc, close all

syms x y xi eta a b x0 y0 E nu t q

PT = [ 1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3 ];
PT_x = diff(PT,x);
PT_y = diff(PT,y);

A = [ ...
      subs(PT,   {x, y}, {-1, -1})
      subs(PT_x, {x, y}, {-1, -1}) % nodo 1
      subs(PT_y, {x, y}, {-1, -1})
      
      subs(PT,   {x, y}, { 1, -1})
      subs(PT_x, {x, y}, { 1, -1}) % nodo 2
      subs(PT_y, {x, y}, { 1, -1})
      
      subs(PT,   {x, y}, { 1,  1})
      subs(PT_x, {x, y}, { 1,  1}) % nodo 3
      subs(PT_y, {x, y}, { 1,  1})
      
      subs(PT,   {x, y}, {-1,  1})
      subs(PT_x, {x, y}, {-1,  1}) % nodo 4
      subs(PT_y, {x, y}, {-1,  1}) ];

N = PT/A;  % N = PT*inv(A);

% reescribo las funciones de forma en terminos de xi y eta
N = subs(N,{x, y},{xi, eta}); 
NN_ = N;

% Reorganizo las funciones de forma en una matriz de 4x3
% donde las filas representan el nodo i
N = [ ...
      N(1)  N(2)  N(3)
      N(4)  N(5)  N(6)       % esto es lo mismo que reshape(N,3,4).'
      N(7)  N(8)  N(9)
      N(10) N(11) N(12) ];

disp('Las funciones de forma son =')
for i = 1:4
   fprintf('N{%d}   = \n',i); pretty(N(i,1));
   fprintf('Nb{%d}  = \n',i); pretty(N(i,2));
   fprintf('Nbb{%d} = \n',i); pretty(N(i,3));
end;

%% Grafico las funciones de forma
XXI  = -1:0.05:1;
EETA = -1:0.05:1;
[XI,ETA] = meshgrid(XXI,EETA);

for i = 1:4
   figure                    % Creo un lienzo
   for j = 1:3
      subplot(1,3,j);        % Divido el lienzo en 3x1 dibujos
      grid on                % creo la rejilla
      hold on;               % Para que no se sobreescriban los graficos

      % con este comando convierto la funcion de forma de tipo simbolico a
      % tipo funcion
      % funcion_texto = char(N(i,j));  % convierto la cadena a texto
      % funcion_texto = strrep(funcion_texto,'*','.*'); % reemplazo * por .*
      % funcion_texto = strrep(funcion_texto,'^','.^'); % reemplazo ^ por .^

      funcion_texto = vectorize(char(N(i,j)));

      NN = inline(funcion_texto,'xi','eta'); % cadena a funcion inline
      % se recomienda aqui mirar la ayuda de la funcion inline
      
      xlabel('\xi', 'FontSize',26); % titulo eje X
      ylabel('\eta','FontSize',26); % titulo eje Y
      switch j % imprimo el titulo
         case 1
            title(sprintf('N_{%d}(\\xi,\\eta)',       i),'FontSize',26);
         case 2
            title(sprintf('(N_b)_{%d}(\\xi,\\eta)',   i),'FontSize',26);
         case 3
            title(sprintf('(N_{bb})_{%d}(\\xi,\\eta)',i),'FontSize',26);
      end
      mesh(XI, ETA, NN(XI,ETA),'LineWidth',2); % malla de alambre
      surf(XI, ETA, NN(XI,ETA));               % superficie
      shading interp         % se interpolan los colores
      alpha 0.3              % opacidad de la superficie
      colormap winter        % mapa de colores a utilizar
      
      axis tight             % ejes apretados
      daspect([1 1 1]);      % similar a axis equal pero en 3D
      view(3);               % vista tridimensional
   end;
end;
drawnow; % vaciar el buffer de graficos de MATLAB antes de continuar

%% Calculo de la matriz de deformacion
BB = cell(4,1);
for i = 1:4
   BB{i} = [...
    -diff(N(i,1),xi,2)/(a^2)           -diff(N(i,2),xi,2)/(a^2)           -diff(N(i,3),xi,2)/(a^2)
    -diff(N(i,1),eta,2)/(b^2)          -diff(N(i,2),eta,2)/(b^2)          -diff(N(i,3),eta,2)/(b^2)      
    -2*diff(diff(N(i,1),xi),eta)/(a*b) -2*diff(diff(N(i,2),xi),eta)/(a*b) -2*diff(diff(N(i,3),xi),eta)/(a*b) ];
end;
B = simple([BB{1} BB{2} BB{3} BB{4}]);

% matriz constitutiva
D = E/(1-nu^2)* [ 1  nu 0
                  nu 1  0
                  0  0  (1-nu)/2 ];
               
Df = (t^3/12)*D; % matriz constitutiva de flexion generalizada

%% Calculo la matriz de rigidez
%{
Recuerde que 
    / dx_dxi   dy_dxi \   / a  0 \
J = |                 | = |      |
    \ dx_deta  dy_dxi /   \ 0  b /

por lo tanto:
%}
det_J = a*b;

tic
disp ('Calculando la matriz de rigidez: espere aproximadamente dos minutos');
K = simple(int(int(B.'*Df*B*det_J, xi, -1, 1), eta, -1, 1))
toc

disp ('Calculando la matriz Df*B');
Df_B = simple(Df*B)

disp ('Calculando el vector de fuerzas nodales equivalentes');
f = int(int(N*q*det_J,xi,-1,1),eta,-1,1);
disp('f = 4*q*a*b*')
disp(simple(f/(4*q*a*b)))
%{
f = 4*q*a*b*
[ 1/4,  1/12,  1/12]
[ 1/4, -1/12,  1/12]
[ 1/4, -1/12, -1/12]
[ 1/4,  1/12, -1/12]
%}

% Calculo de la matriz de masa consistente (FALTA VERIFICAR)
syms rho
disp ('Calculando la matriz de masa consistente');
M = simple(int(int(rho*NN_.'*NN_*det_J,xi,-1,1),eta,-1,1));
disp('M = (a*b*rho)*');
disp(M/(a*b*rho))

return %bye, bye!

%% LA MATRIZ DE RIGIDEZ (copie y pegue en su codigo)
%{
DD = E*t^3/(12*(1-nu^2));   % rigidez a flexion de la placa
Ksimple = simple(K*(a*b)/DD)

Obteniendo:
DD = E*t^3/(12*(1-nu^2));   % rigidez a flexion de la placa
K = DD/(a*b)*[ ...
         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,          (2*nu)/5 + b^2/a^2 + 1/10,          (2*nu)/5 + a^2/b^2 + 1/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             b^2/a^2 - nu/10 + 1/10,      a^2/(2*b^2) - (2*nu)/5 - 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         nu/10 + b^2/(2*a^2) - 1/10,         nu/10 + a^2/(2*b^2) - 1/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      b^2/(2*a^2) - (2*nu)/5 - 1/10,             a^2/b^2 - nu/10 + 1/10
               (2*nu)/5 + b^2/a^2 + 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                 nu,                  nu/10 - b^2/a^2 - 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,              1/10 - b^2/(2*a^2) - nu/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,           b^2/(2*a^2) - (2*nu)/5 - 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0
               (2*nu)/5 + a^2/b^2 + 1/10,                                 nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,           a^2/(2*b^2) - (2*nu)/5 - 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,              1/10 - a^2/(2*b^2) - nu/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,                  nu/10 - a^2/b^2 - 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15
     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             nu/10 - b^2/a^2 - 1/10,      a^2/(2*b^2) - (2*nu)/5 - 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,        - (2*nu)/5 - b^2/a^2 - 1/10,          (2*nu)/5 + a^2/b^2 + 1/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      (2*nu)/5 - b^2/(2*a^2) + 1/10,             a^2/b^2 - nu/10 + 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         1/10 - b^2/(2*a^2) - nu/10,         nu/10 + a^2/(2*b^2) - 1/10
                  b^2/a^2 - nu/10 + 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,             - (2*nu)/5 - b^2/a^2 - 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                -nu,           (2*nu)/5 - b^2/(2*a^2) + 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,              nu/10 + b^2/(2*a^2) - 1/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0
           a^2/(2*b^2) - (2*nu)/5 - 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,               (2*nu)/5 + a^2/b^2 + 1/10,                                -nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,                  nu/10 - a^2/b^2 - 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,              1/10 - a^2/(2*b^2) - nu/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15
 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         1/10 - b^2/(2*a^2) - nu/10,         1/10 - a^2/(2*b^2) - nu/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      (2*nu)/5 - b^2/(2*a^2) + 1/10,             nu/10 - a^2/b^2 - 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,        - (2*nu)/5 - b^2/a^2 - 1/10,        - (2*nu)/5 - a^2/b^2 - 1/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             nu/10 - b^2/a^2 - 1/10,      (2*nu)/5 - a^2/(2*b^2) + 1/10
              nu/10 + b^2/(2*a^2) - 1/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,           (2*nu)/5 - b^2/(2*a^2) + 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,             - (2*nu)/5 - b^2/a^2 - 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                 nu,                  b^2/a^2 - nu/10 + 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0
              nu/10 + a^2/(2*b^2) - 1/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,                  a^2/b^2 - nu/10 + 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,             - (2*nu)/5 - a^2/b^2 - 1/10,                                 nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,           (2*nu)/5 - a^2/(2*b^2) + 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15
     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      b^2/(2*a^2) - (2*nu)/5 - 1/10,             nu/10 - a^2/b^2 - 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         nu/10 + b^2/(2*a^2) - 1/10,         1/10 - a^2/(2*b^2) - nu/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             b^2/a^2 - nu/10 + 1/10,      (2*nu)/5 - a^2/(2*b^2) + 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,          (2*nu)/5 + b^2/a^2 + 1/10,        - (2*nu)/5 - a^2/b^2 - 1/10
           b^2/(2*a^2) - (2*nu)/5 - 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,              1/10 - b^2/(2*a^2) - nu/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,                  nu/10 - b^2/a^2 - 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,               (2*nu)/5 + b^2/a^2 + 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                -nu
                  a^2/b^2 - nu/10 + 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,              nu/10 + a^2/(2*b^2) - 1/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,           (2*nu)/5 - a^2/(2*b^2) + 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,             - (2*nu)/5 - a^2/b^2 - 1/10,                                -nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15 ];

y
Df_B_simple = simple(Df_B*4/DD)
Df_B = DD/4*[...  % = Df*B
  (3*eta*nu*(xi - 1))/b^2 - (3*xi - 3*eta*xi)/a^2,              ((3*xi - 1)*(eta - 1))/a^2,             (nu*(3*eta - 1)*(xi - 1))/b^2,    (3*xi - 3*eta*xi)/a^2 - (3*eta*nu*(xi + 1))/b^2,              ((3*xi + 1)*(eta - 1))/a^2,           -(nu*(3*eta - 1)*(xi + 1))/b^2,  (3*xi + 3*eta*xi)/a^2 + (3*eta*nu*(xi + 1))/b^2,            -((3*xi + 1)*(eta + 1))/a^2,           -(nu*(3*eta + 1)*(xi + 1))/b^2, - (3*xi + 3*eta*xi)/a^2 - (3*eta*nu*(xi - 1))/b^2,            -((3*xi - 1)*(eta + 1))/a^2,             (nu*(3*eta + 1)*(xi - 1))/b^2
 (3*nu*xi*(eta - 1))/a^2 - (3*eta - 3*eta*xi)/b^2,           (nu*(3*xi - 1)*(eta - 1))/a^2,                ((3*eta - 1)*(xi - 1))/b^2, - (3*eta + 3*eta*xi)/b^2 - (3*nu*xi*(eta - 1))/a^2,           (nu*(3*xi + 1)*(eta - 1))/a^2,              -((3*eta - 1)*(xi + 1))/b^2, (3*eta + 3*eta*xi)/b^2 + (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi + 1)*(eta + 1))/a^2,              -((3*eta + 1)*(xi + 1))/b^2,  (3*eta - 3*eta*xi)/b^2 - (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi - 1)*(eta + 1))/a^2,                ((3*eta + 1)*(xi - 1))/b^2
       -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),          ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),       -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b),         ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b) ];

%}
