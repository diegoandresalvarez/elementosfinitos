clear, clc, close all

%% -------------------------------------------------------------------------
%% Funciones de forma del elemento rectangular serendipito de 8 nodos 

% Calculo las funciones de forma unidimensionales
syms xi eta

% Coordenadas de los nodos
%
% Numeracion local:
%      ^ eta
%      |
%      |
%  7---6---5
%  |   |   |
%  8---+---4----> xi
%  |   |   |
%  1---2---3

nod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
    0   -1      %  2
    1   -1      %  3
    1    0      %  4
    1    1      %  5
    0    1      %  6
   -1    1      %  7
   -1    0  ];  %  8

% Se calculan las funciones de forma bidimensionales
xxi  = nod(:,1); 
eeta = nod(:,2);
A = [ ones(8,1) xxi eeta  xxi.^2  xxi.*eeta  eeta.^2  xxi.^2.*eeta  xxi.*eeta.^2 ];

N = cell(8,1);
for i = 1:8
   % se arma el sistema de ecuaciones
   b = zeros(8,1);   b(i) = 1;
   coef_alpha = A\b;
   N{i} = simple([ 1 xi eta xi^2 xi*eta eta^2 xi^2*eta xi*eta^2 ]*coef_alpha);
end

% Imprimo las funciones de forma
fprintf('\nFunciones de forma serendipitas del elemento rectangular de 8 nodos:\n')
for i = 1:8
   fprintf('N%d = %s\n', i, char(N{i}));
end

% Se calculan las derivadas de las funciones de forma con respecto a xi y
% con respecto a eta y se imprimen (para referencias posteriores):
fprintf('\nDerivadas con respecto a xi:\n')
for i = 1:8
   fprintf('dN%d_dxi = %s\n', i, char(simple(diff(N{i},xi))));
end

fprintf('\nDerivadas con respecto a eta:\n')
for i = 1:8
   fprintf('dN%d_deta = %s\n', i, char(simple(diff(N{i},eta))));
end

% Grafico las funciones de forma
XXI  = -1:0.05:1;
EETA = -1:0.05:1;
[XI,ETA] = meshgrid(XXI,EETA);

% Calculo las esferitas
[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i = 1:8
   figure                 % Creo un lienzo
   grid on                % creo la rejilla
   hold on;               % Para que no se sobreescriban los graficos
   
   % con este comando convierto la funcion de forma de tipo simbólico a
   % tipo funcion
%   funcion_texto = char(N{i});  % convierto la cadena a texto  
%   funcion_texto = strrep(funcion_texto,'*','.*'); % reemplazo * por .*
%   funcion_texto = strrep(funcion_texto,'^','.^'); % reemplazo ^ por .^   

   funcion_texto = vectorize(char(N{i}));  % convierto la formula a cadena
   NN = inline(funcion_texto,'xi','eta'); % cadena a funcion inline
   % se recomienda aqui mirar la ayuda de la funcion inline y vectorize
   
   xlabel('\xi','FontSize',26);  % titulo eje X
   ylabel('\eta','FontSize',26); % titulo eje Y
   title(sprintf('N_{%d}(\\xi,\\eta)',i),'FontSize',26); % titulo general
   mesh(XI, ETA, NN(XI,ETA),'LineWidth',2); % malla de alambre
   surf(XI, ETA, NN(XI,ETA));               % superficie
   shading interp         % se interpolan los colores   
   alpha 0.3              % opacidad de la superficie
   colormap winter        % mapa de colores a utilizar
   
   % se grafican las esferitas cada una centrada en (xi_j,eta_j,0)
   for j=1:8
      surf(xsp+nod(j,1), ysp+nod(j,2), zsp+0); 
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1]);      % similar a axis equal pero en 3D
   view(3);               % vista tridimensional
end;

return %bye, bye!
