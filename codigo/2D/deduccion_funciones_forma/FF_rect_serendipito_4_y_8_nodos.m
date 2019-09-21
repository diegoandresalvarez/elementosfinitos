clear, clc, close all

%% Funciones de forma del elemento rectangular serendipito de 4 y 8 nodos 

nno = 8; %  escoja entre {4, 8}.

X = 1; Y = 2;

% coordenadas de los nodos y numeracion local
switch nno
   case 4
      %      ^ eta
      %      |
      %      |
      %  4-------3
      %  |   |   |
      %  ----+---|----> xi
      %  |   |   |
      %  1-------2

      nod = [ ...
      %  xi   eta     % nodo   
         -1   -1      %  1
          1   -1      %  2
          1    1      %  3
         -1    1  ];  %  4
   case 8
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
end      

% se calculan las funciones de forma bidimensionales
xxi  = nod(:, X); 
eeta = nod(:, Y);
switch nno
   case 4
      A = [ ones(4,1) xxi eeta  xxi.*eeta ];
   case 8
      A = [ ones(8,1) xxi eeta  xxi.^2  xxi.*eeta  eeta.^2  xxi.^2.*eeta  xxi.*eeta.^2 ];
end

N = cell(nno,1);
syms xi eta
for i = 1:nno
   % se arma el sistema de ecuaciones
   b = zeros(nno,1);   b(i) = 1;
   coef_alpha = A\b;
   switch nno
      case 4
         N{i} = simplify([ 1 xi eta xi*eta ]*coef_alpha);
      case 8
         N{i} = simplify([ 1 xi eta xi^2 xi*eta eta^2 xi^2*eta xi*eta^2 ]*coef_alpha);
   end
end

%% imprimo las funciones de forma
fprintf('Funciones de forma serendipitas del elemento rectangular de %d nodos:\n', nno)
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

fprintf('\n')
 
%% grafico las funciones de forma
XXI  = linspace(-1, 1, 100);
EETA = linspace(-1, 1, 100);
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