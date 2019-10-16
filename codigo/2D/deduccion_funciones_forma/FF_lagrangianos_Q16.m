clear, clc, close all

%% Funciones de forma del elemento rectangular lagrangiano de 16 nodos 

X = 1; Y = 2;

% Calculo las funciones de forma unidimensionales
syms xi eta
L4_xi = cell(4,1); % contenedor para las funciones de forma (en dir XI)
L4_xi{1} = poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),xi);
L4_xi{2} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),xi);
L4_xi{3} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),xi);
L4_xi{4} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),xi);

L4_eta = cell(4,1); % contenedor para las funciones de forma (en dir ETA)
L4_eta{1} = poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),eta);
L4_eta{2} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),eta);
L4_eta{3} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),eta);
L4_eta{4} = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),eta);

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

nno = size(nod, 1);

% Equivalencia entre coordenada y polinomio
pos = nod;
pos(nod==-1  ) = 1;
pos(nod==-1/3) = 2;
pos(nod== 1/3) = 3;
pos(nod== 1  ) = 4;

% Se calculan las funciones de forma bidimensionales
N = cell(nno,1);
for i = 1:nno
   N{i} = simplify(L4_xi{pos(i,X)}*L4_eta{pos(i,Y)});
end

% Imprimo las funciones de forma
fprintf('\n\nFunciones de forma lagrangianas del elemento rectangular de 16 nodos:\n')
for i = 1:nno
   fprintf('N%d = %s\n', i, char(N{i}))
end

% se calculan las derivadas de las funciones de forma con respecto a xi y
% con respecto a eta y se imprimen (para referencias posteriores):
fprintf('\nDerivadas con respecto a xi:\n')
for i = 1:nno
   fprintf('dN%d_dxi = %s\n',  i, char(simplify(diff(N{i}, xi))))
end

fprintf('\nDerivadas con respecto a eta:\n')
for i = 1:nno
   fprintf('dN%d_deta = %s\n', i, char(simplify(diff(N{i}, eta))))
end

%% Se verifica la condición de cuerpo rígido: sum(N) == 1
suma = 0;
for i = 1:nno
   suma = suma + N{i};
end
fprintf('\nSe verifica la condición de cuerpo rígido: sum(N) == ');
disp(simplify(suma));
 
%% grafico las funciones de forma
XXI  = linspace(-1, 1, 50);
EETA = linspace(-1, 1, 50);
[XI,ETA] = meshgrid(XXI,EETA);

% calculo las esferitas
[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i = 1:nno
   figure                 % creo un lienzo
   grid on                % creo la rejilla
   hold on                % para que no se sobreescriban los graficos
   xlabel('\xi', 'FontSize',16)   % titulo eje X
   ylabel('\eta', 'FontSize',16)  % titulo eje Y
   title(sprintf('N_{%d}(\\xi,\\eta)',i), 'FontSize',20)

   % con este comando convierto la funcion de forma de tipo simbolico a
   % tipo funcion
   NN = matlabFunction(N{i}, 'Vars', [xi, eta]);   
   surf(XI, ETA, NN(XI,ETA))     % superficie
   
   % se grafican las esferitas en cada nodo
   for j=1:nno
      surf(xsp+nod(j,X), ysp+nod(j,Y), zsp+(i==j), 'facecolor', 'k')
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1])       % similar a axis equal pero en 3D
   view(3)                % vista tridimensional
end

%% bye, bye!
return