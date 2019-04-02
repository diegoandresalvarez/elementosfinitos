clear, clc%, close all

%% Programa para calcular los modos de energía nula del sólido del
%% rectángulo serendípito de 8 nodos

X = 1;
Y = 2;
E  = 200;     % modulo de elasticidad del elemento (GPa)
nu = 0.33;    % coeficiente de Poisson
t  = 0.10;    % espesor del elemento (m)

%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
   -((eta + 2*xi)*(eta - 1))/4                  % dN1_dxi
   eta*xi - xi                                  % dN2_dxi
   ((eta - 2*xi)*(eta - 1))/4                   % dN3_dxi
   1/2 - eta^2/2                                % dN4_dxi
   ((eta + 2*xi)*(eta + 1))/4                   % dN5_dxi
   -xi*(eta + 1)                                % dN6_dxi
   -((eta - 2*xi)*(eta + 1))/4                  % dN7_dxi
   eta^2/2 - 1/2                            ];  % dN8_dxi

%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
   -((2*eta + xi)*(xi - 1))/4                   % dN1_deta
   xi^2/2 - 1/2                                 % dN2_deta
   ((xi + 1)*(2*eta - xi))/4                    % dN3_deta
   -eta*(xi + 1)                                % dN4_deta
   ((2*eta + xi)*(xi + 1))/4                    % dN5_deta
   1/2 - xi^2/2                                 % dN6_deta
   -((xi - 1)*(2*eta - xi))/4                   % dN7_deta
   eta*(xi - 1)                             ];  % dN8_deta

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta
n_gl = 3;   % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl]  = gausslegendre_quad(n_gl);

%% matriz constitutiva del elemento para TENSION PLANA
D = [ E/(1-nu^2)     E*nu/(1-nu^2)  0
      E*nu/(1-nu^2)  E/(1-nu^2)     0
      0              0              E/(2*(1+nu)) ];

%% Coordenadas del elemento
xnod = [ ...
%  xi   eta     % nodo
   -1   -1      %  1
    0   -1      %  2
    1   -1      %  3
    1    0      %  4
    1    1      %  5
    0    1      %  6
   -1    1      %  7
   -1    0  ];  %  8

% Calculo las matrices de rigidez y el vector de fuerzas nodales
% equivalentes del elemento
K = zeros(16);
B = cell(n_gl,n_gl); % contenedor para las matrices de deformacion
det_J = zeros(n_gl,n_gl); % en esta matriz se almacenaran los Jacobianos
for p = 1:n_gl
   for q = 1:n_gl
      xi_gl  = x_gl(p);
      eta_gl = x_gl(q);
           
      % Se evaluan las derivadas de las funciones de forma en los puntos
      % de integracion de Gauss-Legendre
      ddN_dxi  = dN_dxi (xi_gl, eta_gl);       xe = xnod(:,X);
      ddN_deta = dN_deta(xi_gl, eta_gl);       ye = xnod(:,Y);
      
      dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
      dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
      
      % Se ensambla la matriz Jacobiana del elemento
      J = [ dx_dxi   dy_dxi
            dx_deta  dy_deta ];
      
      % Se calcula el determinante del Jacobiano
      det_J(p,q) = det(J);
      
      B{p,q} = zeros(3,2*8);
      for i = 1:8        
         % Se ensambla la matriz de deformacion del elemento B
         dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_J(p,q);
         dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_J(p,q);
         B{p,q}(:,[2*i-1 2*i]) = [ dNi_dx 0          % aqui se ensambla
                                   0      dNi_dy     % y asigna la matriz
                                   dNi_dy dNi_dx ];  % B_i
      end;
      
      % se arma la matriz de rigidez del elemento e
      K = K + B{p,q}'*D*B{p,q}*det_J(p,q)*t*w_gl(p)*w_gl(q);     
   end;
end;

%% Se calculan los valores y vectores propios de la matriz K
[evec,eval] = eig(K);
[eval,idx]  = sort(diag(eval));
evec = evec(:,idx);

%% Se imprimen los vectores propios (recuerde que los modos de energía nula
%% son aquellos para los cuales los valores propios son cero
figure
modo = cell(16,1);
for i = 1:16
   modo{i} = reshape(evec(:,i),2,8)' + xnod;
   subplot(4,4,i)
   plot(xnod([1:8 1],X),xnod([1:8 1],Y),'b');
   hold on
   plot(modo{i}([1:8 1],X), modo{i}([1:8 1],Y),'r');
   axis equal
   axis([-2 2 -2 2]);
   title(sprintf('\\lambda_%d = %d',i,eval(i)),'FontSize',16);
end

ax = axes('Position', [0,0,1,1], 'Visible', 'off');
tx = text(0.4, 0.95, sprintf('Puntos de integración = %d', n_gl));
set(tx, 'FontWeight', 'Bold', 'FontSize', 20);

