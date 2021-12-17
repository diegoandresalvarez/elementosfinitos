%% Funciones de forma de los EFs triangulares de 3, 6 y 10 nodos 

clear, clc, close all

syms xi L1 L2 L3

%% Funciones de forma del EF triangular de 3 nodos
disp('Funciones de forma de triangulos de 3 nodos')

T3.coord = [1 0 0
            0 1 0
            0 0 1];

IJK = round(1*T3.coord);
%       [1 0 0
%        0 1 0
%        0 0 1];
    
T3.N = cell(3,1);
for i = 1:3
   switch IJK(i,1)
      case 0, lI = 1;
      case 1, lI = calc_N([0 1], [0 1], L1);
   end
   switch IJK(i,2)
      case 0, lJ = 1;
      case 1, lJ = calc_N([0 1], [0 1], L2);
   end
   switch IJK(i,3)
      case 0, lK = 1;
      case 1, lK = calc_N([0 1], [0 1], L3);
   end
     
   T3.N{i} = lI*lJ*lK; % = lI^i(L1) * lJ^i(L2) * lK^i(L3)
   fprintf('N{%d} = ',i); disp(T3.N{i});
end
%EF = T3;

%% Funciones de forma del EF triangular de 6 nodos
fprintf('\nFunciones de forma de triangulos de 6 nodos\n')

T6.coord = [1    0    0
            0    1    0
            0    0    1
            1/2  1/2  0
            0    1/2  1/2
            1/2  0    1/2];

IJK = round(2*T6.coord);
%       [2 0 0
%        0 2 0
%        0 0 2
%        1 1 0
%        0 1 1
%        1 0 1];
    
T6.N = cell(6,1);
for i = 1:6  
   switch IJK(i,1)
      case 0, lI = 1;
      case 1, lI = calc_N([0  1/2   ], [0 1  ], L1);
      case 2, lI = calc_N([0  1/2  1], [0 0 1], L1);
   end
   switch IJK(i,2)
      case 0, lJ = 1;
      case 1, lJ = calc_N([0  1/2   ], [0 1  ], L2);
      case 2, lJ = calc_N([0  1/2  1], [0 0 1], L2);
   end
   switch IJK(i,3)
      case 0, lK = 1;
      case 1, lK = calc_N([0  1/2   ], [0 1  ], L3);
      case 2, lK = calc_N([0  1/2  1], [0 0 1], L3);
   end
 
   T6.N{i} = simplify(lI*lJ*lK); % = lI^i(L1) * lJ^i(L2) * lK^i(L3)
   fprintf('N{%d} = ',i); disp(T6.N{i});
end
%EF = T6;

%% Funciones de forma del EF triangular de 10 nodos
fprintf('\nFunciones de forma de triangulos de 10 nodos\n')

T10.coord = [1    0    0
             0    1    0
             0    0    1
             2/3  1/3  0
             1/3  2/3  0
             0    2/3  1/3
             0    1/3  2/3
             1/3  0    2/3
             2/3  0    1/3             
             1/3  1/3  1/3];

IJK = round(3*T10.coord);
%       [3 0 0
%        0 3 0
%        0 0 3
%        2 1 0
%        1 2 0
%        0 2 1
%        0 1 2
%        1 0 2
%        2 0 1
%        1 1 1];
    
T10.N = cell(10,1);
for i = 1:10  
   switch IJK(i,1)
      case 0, lI = 1;
      case 1, lI = calc_N([0  1/3        ], [0 1    ], L1);
      case 2, lI = calc_N([0  1/3  2/3   ], [0 0 1  ], L1);
      case 3, lI = calc_N([0  1/3  2/3  1], [0 0 0 1], L1);
   end
   switch IJK(i,2)
      case 0, lJ = 1;
      case 1, lJ = calc_N([0  1/3        ], [0 1    ], L2);
      case 2, lJ = calc_N([0  1/3  2/3   ], [0 0 1  ], L2);
      case 3, lJ = calc_N([0  1/3  2/3  1], [0 0 0 1], L2);
   end
   switch IJK(i,3)
      case 0, lK = 1;
      case 1, lK = calc_N([0  1/3        ], [0 1    ], L3);
      case 2, lK = calc_N([0  1/3  2/3   ], [0 0 1  ], L3);
      case 3, lK = calc_N([0  1/3  2/3  1], [0 0 0 1], L3);
   end
     
   T10.N{i} = simplify(lI*lJ*lK); % = lI^i(L1) * lJ^i(L2) * lK^i(L3)
   fprintf('N{%d} = ',i); disp(T10.N{i});
end
EF = T10;

%% Se grafican las funciones de forma
EF.nno = size(EF.coord, 1);

LL2 = 0:0.05:1;
LL3 = 0:0.05:1;
[L2,L3] = meshgrid(LL2,LL3);
L1 = 1 - L2 - L3;
L1 = round(100*L1)/100;   % redondear bien!
L1(L1<0) = NaN;

% coordenadas del triangulo
x = [0 1.0 0.5]; y = [0 0 sqrt(0.75)];
X = L1*x(1) + L2*x(2) + L3*x(3);
Y = L1*y(1) + L2*y(2) + L3*y(3);
X = X(:); 
Y = Y(:); 
isnanX = isnan(X);
X(isnanX) = [];
Y(isnanX) = [];

TRI = delaunay(X,Y);
numtriang = size(TRI, 1);   % numero de triangulos
% Se eliminan de la lista de triangulos aquellos cuya area sea negativa o 
% igual a cero
Area = zeros(numtriang,1);
for i = 1:numtriang
   Area(i) = 0.5*det([ 1 X(TRI(i,1)) Y(TRI(i,1))      %Area del EF e
                       1 X(TRI(i,2)) Y(TRI(i,2))
                       1 X(TRI(i,3)) Y(TRI(i,3))]);
end
TRI(Area < 1e-3,:) = [];

[xsp,ysp,zsp] = sphere;
xsp = 0.025*xsp;
ysp = 0.025*ysp;
zsp = 0.025*zsp;

for i=1:EF.nno
   % creo una funcion que pueda evaluar numericamente
   N = matlabFunction(EF.N{i}, 'Vars', [sym('L1'),sym('L2'),sym('L3')]);

   Z = N(L1,L2,L3);
   Z = Z(:);
   Z = Z(~isnanX);
   figure
   trimesh(TRI,X,Y,Z,'LineWidth',2);
   hold on
   trisurf(TRI,X,Y,Z);   
   alpha 0.3
   shading interp
   colormap winter
%   axis([-0.1 1.1 -0.1 0.96 -0.4 1.1])
   axis tight
   for j=1:EF.nno
      cxsp = EF.coord(j,:)*x';
      cysp = EF.coord(j,:)*y';
      surf(xsp+cxsp, ysp+cysp, zsp+0);  % sphere centered at (x(1),y(1),0)
   end
   daspect([1 1 1]);
   title(sprintf('N_{%d} = %s',i,char(EF.N{i})),'FontSize',20);
%   print('-dpdf',sprintf('%d.pdf',i));
end

%% Se verifica la condicion de cuerpo rigido: sum(N) == 1
syms a b
suma = subs(sum([EF.N{:}]), {'L1', 'L2', 'L3'}, {1-a-b, a, b});
fprintf('\nSe verifica la condicion de cuerpo rigido: sum(N) == ');
disp(simplify(suma));

%% Bye, bye!
return;

%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamanios de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamanio');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximacion numerica, haciendo los
    % coeficientes demasiado pequenios igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end
