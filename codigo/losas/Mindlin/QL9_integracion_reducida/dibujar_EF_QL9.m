function dibujar_EF_QL9(xe, ye, N, ae, te, esc_w, esc_t)
%% Dibuja un EF de losa QL9 deformado
% xe, ye  nodos del borde de la losa
% N       funciones de forma de la losa
% ae      desplazamientos nodales en los nodos del borde de la losa
% te      espesor de la losa
% esc_w   escalamiento del desplazamiento vertical
% esc_t   escalamiento de los giros

%%
%{
% EJEMPLO PARA EL ELEMENTO RECTANGULAR ISOPARAMETRICO SERENDIPITO DE 8
% NODOS:
% Numeracion local:
%      ^ eta
%      |
%      |
%  7---6---5
%  |   |   |
%  8---+---4----> xi
%  |   |   |
%  1---2---3

xnod = [ ...
%   x     y       % nodo   
   -1.5  -1       %  1
    0    -0.9     %  2
    1    -1       %  3
    1.3   0       %  4
    1     1.2     %  5
    0.2   1.1     %  6
   -1     0.7     %  7
   -1.3   0   ];  %  8

xe = xnod(:,1); ye = xnod(:,2);

Nforma = @(xi,eta) [ ...
-((eta - 1)*(xi - 1)*(eta + xi + 1))/4       % N1
((xi^2 - 1)*(eta - 1))/2                     % N2
((eta - 1)*(xi + 1)*(eta - xi + 1))/4        % N3
-((eta^2 - 1)*(xi + 1))/2                    % N4
((eta + 1)*(xi + 1)*(eta + xi - 1))/4        % N5
-((xi^2 - 1)*(eta + 1))/2                    % N6
((eta + 1)*(xi - 1)*(xi - eta + 1))/4        % N7
((eta^2 - 1)*(xi - 1))/2                 ];  % N8

ae = 0.1*randn(24,1);
%ae(2:3:end) = 0.1*ae(2:3:end);
%ae(3:3:end) = 0.1*ae(3:3:end);

te = 0.2;

figure
grid on
hold on
view(3)
esc_w = 3; esc_t = 3;
dibujar_EF_QL9(xe, ye, Nforma, ae, te, esc_w, esc_t)
daspect([1 1 2]);
axis tight

%}

%% Se definen algunas constantes
w_ = 1;  tx_ = 2;  ty_ = 3; % lectura del codigo

nno = size(xe,1);           % numero de nodos del elemento finito

%% Se calcula la geometria isoparametrica
xi  = [               -1:0.1:1   +ones(1,length(-1:0.1:1))                1:-0.1:-1   -ones(1,length(1:-0.1:-1)) ];
eta = [-ones(1,length(-1:0.1:1))               -1:0.1:1    +ones(1,length(1:-0.1:-1))                1:-0.1:-1   ];
n = length(xi);  % number of points of vector xi

x = zeros(n,1);
y = zeros(n,1);
for j = 1:n
   x(j) = sum(N(xi(j), eta(j)) .* xe);
   y(j) = sum(N(xi(j), eta(j)) .* ye);
end

%% Se grafica el elemento original
plot3(x, y, repmat(-te/2,size(x)), 'b-'); % cara inferior
plot3(x, y, zeros(       size(x)), 'b-'); % plano medio
plot3(x, y, repmat( te/2,size(x)), 'b-'); % cara superior

%% Se grafican las lineas verticales originales
for j = 1:nno
   plot3([ xe(j) xe(j) ], [ ye(j) ye(j) ], [ -te/2 te/2 ], 'b*-');
end

%% Se calcula el elemento deformado
mov = reshape(ae,3,nno)';  % se reorganiza el vector de movimientos nodales
mov(:,w_)        = esc_w*mov(:,w_);        % se escala la deform vertical
mov(:,[tx_ ty_]) = esc_t*mov(:,[tx_ ty_]); % se escalan los angulos

w  = zeros(n,1);
tx = zeros(n,1);
ty = zeros(n,1);
for j = 1:n
   w(j)  = sum(N(xi(j), eta(j)) .* mov(:,w_));
   tx(j) = sum(N(xi(j), eta(j)) .* mov(:,tx_));
   ty(j) = sum(N(xi(j), eta(j)) .* mov(:,ty_));   
end
z =  te/2; u_sup = -z*tx; v_sup = -z*ty; w_sup = w+z;
z =     0; u_mid = -z*tx; v_mid = -z*ty; w_mid = w+z;
z = -te/2; u_inf = -z*tx; v_inf = -z*ty; w_inf = w+z;

%% Se grafica el elemento deformado
plot3(x+u_sup, y+v_sup, w_sup, 'r'); % cara superior
plot3(x+u_mid, y+v_mid, w_mid, 'r'); % plano medio
plot3(x+u_inf, y+v_inf, w_inf, 'r'); % cara inferior

%% Se grafican las lineas verticales deformadas
z =  te/2; u_sup = -z*mov(:,tx_); v_sup = -z*mov(:,ty_); w_sup = mov(:,w_)+z;
z = -te/2; u_inf = -z*mov(:,tx_); v_inf = -z*mov(:,ty_); w_inf = mov(:,w_)+z;

for j = 1:nno
   plot3([ xe(j)+u_inf(j) xe(j)+u_sup(j) ], [ ye(j)+v_inf(j) ye(j)+v_sup(j) ], [ w_inf(j) w_sup(j) ], 'r*-');
end

%%
[xi, eta] = meshgrid(-1:0.1:1);
w = zeros(21);
x = zeros(21);
y = zeros(21);
for i = 1:21
   for j = 1:21
      w(i,j) = sum(N(xi(i,j), eta(i,j)) .* mov(:,w_));
      x(i,j) = sum(N(xi(i,j), eta(i,j)) .* xe);
      y(i,j) = sum(N(xi(i,j), eta(i,j)) .* ye);
   end
end

surf(x,y,w);

%%
return, % bye, bye !!
