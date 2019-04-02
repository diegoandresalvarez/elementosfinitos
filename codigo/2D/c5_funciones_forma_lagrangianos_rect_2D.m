clear, clc, close all

%% -------------------------------------------------------------------------
%% Funciones de forma del elemento rectangular lagrangiano de 16 nodos 

% Calculo las funciones de forma unidimensionales
syms xi eta
L4_xi = cell(4,1); % contenedor para las funciones de forma (en dir XI)
L4_xi{1} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),xi));
L4_xi{2} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),xi));
L4_xi{3} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),xi));
L4_xi{4} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),xi));

L4_eta = cell(4,1); % contenedor para las funciones de forma (en dir ETA)
L4_eta{1} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),eta));
L4_eta{2} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),eta));
L4_eta{3} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),eta));
L4_eta{4} = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),eta));

% Coordenadas de los nodos
%
% Numeración local:
%        ^ eta
%        |
%        |
% 10---9---8---7
%  |   |   |   |
% 11--16--15---6
%  |   |   |   |----> xi
% 12--13--14---5
%  |   |   |   |
%  1---2---3---4

nod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
   -1/3 -1      %  2
    1/3 -1      %  3
    1   -1      %  4
    1   -1/3    %  5
    1    1/3    %  6
    1    1      %  7
    1/3  1      %  8
   -1/3  1      %  9
   -1    1      % 10
   -1    1/3    % 11
   -1   -1/3    % 12
   -1/3 -1/3    % 13
    1/3 -1/3    % 14
    1/3  1/3    % 15
   -1/3  1/3 ]; % 16

% Equivalencia entre coordenada y polinomio
pos = nod;
pos(nod==-1  ) = 1;
pos(nod==-1/3) = 2;
pos(nod== 1/3) = 3;
pos(nod== 1  ) = 4;

% Se calculan las funciones de forma bidimensionales
N = cell(16,1);
for i = 1:16
   N{i} = simple(L4_xi{pos(i,1)}*L4_eta{pos(i,2)});
end

% Imprimo las funciones de forma
fprintf('\n\nFunciones de forma lagrangianas del elemento rectangular de 16 nodos:\n')
for i = 1:16
   fprintf('\n\nN%d = ',i), pretty(N{i});
end

% Se calculan las derivadas de las funciones de forma con respecto a xi y
% con respecto a eta y se imprimen (para referencias posteriores):
fprintf('\n\nDerivadas con respecto a xi:\n')
for i = 1:16
   fprintf('\n\ndN%d_dxi = ',i), pretty(simple(diff(N{i},xi)));
end

fprintf('\n\nDerivadas con respecto a eta:\n')
for i = 1:16
   fprintf('\n\ndN%d_deta = ',i), pretty(simple(diff(N{i},eta)));
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

for i = 1:16
   figure                 % Creo un lienzo
   grid on                % creo la rejilla
   hold on;               % Para que no se sobreescriban los graficos
   
   % con este comando convierto la función de forma de tipo simbólico a
   % tipo función
   funcion_texto = char(N{i});  % convierto la cadena a texto
   funcion_texto = strrep(funcion_texto,'*','.*'); % reemplazo * por .*
   funcion_texto = strrep(funcion_texto,'^','.^'); % reemplazo ^ por .^   
   NN = inline(funcion_texto,'xi','eta'); % cadena a funcion inline
   % se recomienda aqui mirar la ayuda de la funcion inline
   
   xlabel('\xi','FontSize',26);  % titulo eje X
   ylabel('\eta','FontSize',26); % titulo eje Y
   title(sprintf('N_{%d}(\\xi,\\eta)',i),'FontSize',26); % titulo general
   mesh(XI, ETA, NN(XI,ETA),'LineWidth',2); % malla de alambre
   surf(XI, ETA, NN(XI,ETA));               % superficie
   shading interp         % se interpolan los colores   
   alpha 0.3              % opacidad de la superficie
   colormap winter        % mapa de colores a utilizar
   
   % se grafican las esferitas cada una centrada en (xi_j,eta_j,0)
   for j=1:16
      surf(xsp+nod(j,1), ysp+nod(j,2), zsp+0); 
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1]);      % similar a axis equal pero en 3D
   view(3);               % vista tridimensional
end;

return %bye, bye!
