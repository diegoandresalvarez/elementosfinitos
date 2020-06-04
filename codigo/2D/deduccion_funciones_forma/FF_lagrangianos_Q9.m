clear, clc, close all

%% Funciones de forma del elemento rectangular lagrangiano de 9 nodos 

X = 1; Y = 2;

% Calculo las funciones de forma unidimensionales
syms xi eta
L3_xi = cell(3,1); % contenedor para las funciones de forma (en dir XI)
L3_xi{1} = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);
L3_xi{2} = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);
L3_xi{3} = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);

L3_eta = cell(3,1); % contenedor para las funciones de forma (en dir ETA)
L3_eta{1} = poly2sym(polyfit([-1 0 1],[1 0 0],2),eta);
L3_eta{2} = poly2sym(polyfit([-1 0 1],[0 1 0],2),eta);
L3_eta{3} = poly2sym(polyfit([-1 0 1],[0 0 1],2),eta);

% Coordenadas de los nodos
%
% Numeracion local:
%     ^ eta
%     |
%     |
% 7---6---5
% |   |   |
% 8---9---4----> xi
% |   |   |
% 1---2---3

nod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
    0   -1      %  2
    1   -1      %  3
    1    0      %  4
    1    1      %  5
    0    1      %  6
   -1    1      %  7
   -1    0      %  8
    0    0 ];   %  9

nno = size(nod, 1);

% Equivalencia entre coordenada y polinomio
pos = nod;
pos(nod==-1) = 1;
pos(nod== 0) = 2;
pos(nod== 1) = 3;

% Se calculan las funciones de forma bidimensionales
N = cell(nno,1);
for i = 1:nno
   N{i} = simplify(L3_xi{pos(i,1)}*L3_eta{pos(i,2)});
end

% Imprimo las funciones de forma
fprintf('Funciones de forma lagrangianas del elemento rectangular de 9 nodos:\n')
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

%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
suma = 0;
for i = 1:nno
   suma = suma + N{i};
end
fprintf('\nSe verifica la condicion de cuerpo rigido: sum(N) == ');
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
