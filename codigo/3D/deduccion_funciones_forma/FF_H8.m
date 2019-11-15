clear, clc, close all

%% -------------------------------------------------------------------------
%% Funciones de forma del elemento hexahedrico lagrangiano de 8 nodos

% Coordenadas de los nodos
%
% Numeracion local:
%      ^ eta             para zeta = -1
%      |
%      |
%  4-------3
%  |   |   |
%  |---+---|----> xi
%  |   |   |
%  1-------2
%
%
% Numeracion local:
%      ^ eta             para zeta = +1
%      |
%      |
%  8-------7
%  |   |   |
%  |---+---|----> xi
%  |   |   |
%  5-------6

nod = [ ...
%  xi   eta  zeta   nodo
   -1   -1   -1    % 1
    1   -1   -1    % 2
    1    1   -1    % 3
   -1    1   -1    % 4
   -1   -1    1    % 5
    1   -1    1    % 6
    1    1    1    % 7
   -1    1    1 ]; % 8

%% Calculo las funciones de forma unidimensionales
syms xi eta zeta
L2_xi = cell(2,1); % contenedor para las funciones de forma (en dir XI)
L2_xi{1} = poly2sym(polyfit([-1 1],[1 0],1), xi);
L2_xi{2} = poly2sym(polyfit([-1 1],[0 1],1), xi);

L2_eta = cell(2,1); % contenedor para las funciones de forma (en dir ETA)
L2_eta{1} = poly2sym(polyfit([-1 1],[1 0],1), eta);
L2_eta{2} = poly2sym(polyfit([-1 1],[0 1],1), eta);

L2_zeta = cell(2,1); % contenedor para las funciones de forma (en dir ZETA)
L2_zeta{1} = poly2sym(polyfit([-1 1],[1 0],1), zeta);
L2_zeta{2} = poly2sym(polyfit([-1 1],[0 1],1), zeta);

%% Equivalencia entre coordenada y polinomio
pos = nod;
pos(nod==-1) = 1;
pos(nod== 1) = 2;

%% Se calculan las funciones de forma tridimensionales
N = cell(8,1);
for i = 1:8
   N{i} = simplify(L2_xi{pos(i,1)} * L2_eta{pos(i,2)} * L2_zeta{pos(i,3)});
end

ev = zeros(8,1);
for i = 1:8
   NN = matlabFunction(N{i}, 'Vars', {'xi','eta','zeta'});
   for j = 1:8
      ev(j) = NN(nod(j,1), nod(j,2), nod(j,3));
   end
   fprintf('%d = \n', i)
   disp(ev)
end

%% Imprimo las funciones de forma
fprintf('\n\nFunciones de forma lagrangianas del elemento hexahedrico de 8 nodos:\n')
for i = 1:8
   fprintf('\nN{%d} = %s', i, N{i});
end

%% Se calculan las derivadas de las funciones de forma con respecto a xi y
% con respecto a eta y se imprimen (para referencias posteriores):
fprintf('\n\nDerivadas con respecto a xi:\n')
for i = 1:8
   fprintf('\ndN%d_dxi = %s',i, simplify(diff(N{i},xi)));
end

fprintf('\n\nDerivadas con respecto a eta:\n')
for i = 1:8
   fprintf('\ndN%d_deta = %s',i, simplify(diff(N{i},eta)));
end

fprintf('\n\nDerivadas con respecto a zeta:\n')
for i = 1:8
   fprintf('\ndN%d_dzeta = %s',i, simplify(diff(N{i},zeta)));
end

%% Se verifica la condición de cuerpo rígido: sum(N) == 1
suma = 0;
for i = 1:8
   suma = suma + N{i};
end
fprintf('\nSe verifica la condición de cuerpo rígido: sum(N) == ');
disp(simplify(suma));

%% Grafico las funciones de forma
XXI   = -1:0.1:1;
EETA  = -1:0.1:1;
ZZETA = -1:0.1:1;
[XI,ETA,ZETA] = meshgrid(XXI,EETA,ZZETA);

% Calculo las esferitas
[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i = 1:8
   figure                 % Creo un lienzo
   grid on                % creo la rejilla
   hold on;               % Para que no se sobreescriban los graficos 
   xlabel('\xi',  'FontSize',20); % titulo eje X
   ylabel('\eta', 'FontSize',20); % titulo eje Y
   zlabel('\zeta','FontSize',20); % titulo eje Z
   title(sprintf('N_{%d}(\\xi,\\eta,\\zeta)',i),'FontSize',20); % titulo general

   NN = matlabFunction(N{i}, 'Vars', {'xi','eta','zeta'});   
   xslice = [-1 0 1]; yslice = [-1 0 1]; zslice = [-1 0 1];
   slice(XI,ETA,ZETA, NN(XI,ETA,ZETA), xslice,yslice,zslice);

   shading interp         % se interpolan los colores   
   alpha 0.3              % opacidad de la superficie
   
   % se grafican las esferitas cada una centrada en (xi_j,eta_j,zeta_j)
   for j=1:8
      surf(xsp+nod(j,1), ysp+nod(j,2), zsp+nod(j,3), 'facecolor', 'k')
      h = text(nod(j,1), nod(j,2), nod(j,3), num2str(j));
      set(h,'Color', [1 0 0], 'FontSize',16);
   end 
   
   axis tight             % ejes apretados
   daspect([1 1 1]);      % similar a axis equal pero en 3D
   view(3);               % vista tridimensional
   colorbar
end

%% Calculo de la matriz de forma
NN = sym(zeros(3,3*8));
for i = 1:8
	% Se ensambla la matriz de funciones de forma N
   NN(:,3*i-[2 1 0]) = diag([N{i} N{i} N{i}]);
end

%% Calculo de la matriz de masa M
syms rho V
det_J = V;   % el determinante del Jacobiano es el volumen del hexahedro
M = simplify(int(int(int(rho*(NN.'*NN)*det_J, xi,-1,1), eta,-1,1), zeta,-1,1));
fprintf('\n\nM = rho*V/27*\n'); pretty(27*M/(rho*V));

return %bye, bye!
