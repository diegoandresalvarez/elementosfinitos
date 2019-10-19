clear,clc,close all

xnod = zeros(16,1);
ynod = zeros(16,1);
figure(1)
title('De click con el mouse 16 veces')
axis([-5 5 -5 5]);
hold on;
for i = 1:16
   [xnod(i),ynod(i)] = ginput(1);
   plot(xnod(i), ynod(i), 'ko','MarkerSize',12);  
   plot(xnod(i), ynod(i), 'kx','MarkerSize',12);
end
close(1)

delta = 0.05;
xxi  = -1:delta:1;
eeta = -1:delta:1;
n = length(xxi);
[xi,eta] = meshgrid(xxi,eeta);

figure
subplot(2,2,1);
hold on;
for i = 1:n
   h1 = plot(xi(:,i),eta(:,i));
   h2 = plot(xi(i,:),eta(i,:));
   if i==1 || i==n
      set(h1, 'LineWidth', 4);
      set(h2, 'LineWidth', 4);
   end
end
axis equal; axis([-1.1 1.1 -1.1 1.1])

% Funciones de forma del elemento lagrangiano plano de 16 nodos (cuadratico)
%
% Numeracion local:
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

Ni1 = (1/16).*(xi-1).*(1-9.*xi.^2);
Ni2 = (9/16).*(1-xi.^2).*(1-3.*xi);
Ni3 = (9/16).*(1-xi.^2).*(1+3.*xi);
Ni4 = (1/16).*(xi+1).*(9.*xi.^2-1);

Nj1 = (1/16).*(eta-1).*(1-9.*eta.^2);
Nj2 = (9/16).*(1-eta.^2).*(1-3.*eta);
Nj3 = (9/16).*(1-eta.^2).*(1+3.*eta);
Nj4 = (1/16).*(eta+1).*(9.*eta.^2-1);

N = cell(16,1);
N{10} = Ni1.*Nj4;   N{9}  = Ni2.*Nj4;   N{8}  = Ni3.*Nj4;   N{7}  = Ni4.*Nj4;
N{11} = Ni1.*Nj3;   N{16} = Ni2.*Nj3;   N{15} = Ni3.*Nj3;   N{6}  = Ni4.*Nj3;
N{12} = Ni1.*Nj2;   N{13} = Ni2.*Nj2;   N{14} = Ni3.*Nj2;   N{5}  = Ni4.*Nj2;
N{1}  = Ni1.*Nj1;   N{2}  = Ni2.*Nj1;   N{3}  = Ni3.*Nj1;   N{4}  = Ni4.*Nj1;

x = zeros(n);
y = zeros(n);
for i = 1:16
   x = x + N{i}*xnod(i);
   y = y + N{i}*ynod(i);
end
xinod  = [-1 -1/3 1/3 1   1   1   1 1/3 -1/3 -1 -1  -1   -1/3  1/3 1/3 -1/3];
etanod = [-1 -1   -1  -1 -1/3 1/3 1 1    1    1 1/3 -1/3 -1/3 -1/3 1/3  1/3];
for i = 1:16
   plot(xinod(i), etanod(i), 'r*','MarkerSize',12, 'LineWidth', 4);  
end


subplot(2,2,[2 4]);
hold on;
for i = 1:n
   h1 = plot(x(:,i),y(:,i));
   h2 = plot(x(i,:),y(i,:));
   if i==1 || i==n
      set(h1, 'LineWidth', 4);
      set(h2, 'LineWidth', 4);
   end
end
for i = 1:16
   plot(xnod(i), ynod(i), 'r*','MarkerSize',12, 'LineWidth', 4);  
end
axis equal tight

dN_dxi{1} = ((- 27.*xi.^2 + 18.*xi + 1).*(- 9.*eta.^3 + 9.*eta.^2 + eta - 1))/256;
dN_dxi{2} = -(9.*(- 9.*xi.^2 + 2.*xi + 3).*(- 9.*eta.^3 + 9.*eta.^2 + eta - 1))/256;
dN_dxi{3} = -(9.*(9.*xi.^2 + 2.*xi - 3).*(- 9.*eta.^3 + 9.*eta.^2 + eta - 1))/256;
dN_dxi{4} = ((27.*xi.^2 + 18.*xi - 1).*(- 9.*eta.^3 + 9.*eta.^2 + eta - 1))/256;
dN_dxi{5} = (9.*(3.*eta - 1).*(eta - 1).*(eta + 1).*(27.*xi.^2 + 18.*xi - 1))/256;
dN_dxi{6} = (9.*(27.*xi.^2 + 18.*xi - 1).*(- 3.*eta.^3 - eta.^2 + 3.*eta + 1))/256;
dN_dxi{7} = -((27.*xi.^2 + 18.*xi - 1).*(- 9.*eta.^3 - 9.*eta.^2 + eta + 1))/256;
dN_dxi{8} = (9.*(9.*xi.^2 + 2.*xi - 3).*(- 9.*eta.^3 - 9.*eta.^2 + eta + 1))/256;
dN_dxi{9} = (9.*(- 9.*xi.^2 + 2.*xi + 3).*(- 9.*eta.^3 - 9.*eta.^2 + eta + 1))/256;
dN_dxi{10} = -((- 27.*xi.^2 + 18.*xi + 1).*(- 9.*eta.^3 - 9.*eta.^2 + eta + 1))/256;
dN_dxi{11} = (9.*(eta - 1).*(3.*eta + 1).*(eta + 1).*(27.*xi.^2 - 18.*xi - 1))/256;
dN_dxi{12} = -(9.*(eta - 1).*(3.*eta - 1).*(eta + 1).*(27.*xi.^2 - 18.*xi - 1))/256;
dN_dxi{13} = (81.*(eta - 1).*(3.*eta - 1).*(eta + 1).*(9.*xi.^2 - 2.*xi - 3))/256;
dN_dxi{14} = (81.*(9.*xi.^2 + 2.*xi - 3).*(- 3.*eta.^3 + eta.^2 + 3.*eta - 1))/256;
dN_dxi{15} = (81.*(3.*eta + 1).*(eta - 1).*(eta + 1).*(9.*xi.^2 + 2.*xi - 3))/256;
dN_dxi{16} = -(81.*(eta - 1).*(3.*eta + 1).*(eta + 1).*(9.*xi.^2 - 2.*xi - 3))/256;
dN_deta{1} = ((- 27.*eta.^2 + 18.*eta + 1).*(- 9.*xi.^3 + 9.*xi.^2 + xi - 1))/256;
dN_deta{2} = -(9.*(xi - 1).*(3.*xi - 1).*(xi + 1).*(27.*eta.^2 - 18.*eta - 1))/256;
dN_deta{3} = (9.*(xi - 1).*(3.*xi + 1).*(xi + 1).*(27.*eta.^2 - 18.*eta - 1))/256;
dN_deta{4} = -((- 27.*eta.^2 + 18.*eta + 1).*(- 9.*xi.^3 - 9.*xi.^2 + xi + 1))/256;
dN_deta{5} = (9.*(- 9.*eta.^2 + 2.*eta + 3).*(- 9.*xi.^3 - 9.*xi.^2 + xi + 1))/256;
dN_deta{6} = (9.*(9.*eta.^2 + 2.*eta - 3).*(- 9.*xi.^3 - 9.*xi.^2 + xi + 1))/256;
dN_deta{7} = -((27.*eta.^2 + 18.*eta - 1).*(- 9.*xi.^3 - 9.*xi.^2 + xi + 1))/256;
dN_deta{8} = (9.*(27.*eta.^2 + 18.*eta - 1).*(- 3.*xi.^3 - xi.^2 + 3.*xi + 1))/256;
dN_deta{9} = (9.*(3.*xi - 1).*(xi - 1).*(xi + 1).*(27.*eta.^2 + 18.*eta - 1))/256;
dN_deta{10} = ((27.*eta.^2 + 18.*eta - 1).*(- 9.*xi.^3 + 9.*xi.^2 + xi - 1))/256;
dN_deta{11} = -(9.*(9.*eta.^2 + 2.*eta - 3).*(- 9.*xi.^3 + 9.*xi.^2 + xi - 1))/256;
dN_deta{12} = -(9.*(- 9.*eta.^2 + 2.*eta + 3).*(- 9.*xi.^3 + 9.*xi.^2 + xi - 1))/256;
dN_deta{13} = (81.*(xi - 1).*(3.*xi - 1).*(xi + 1).*(9.*eta.^2 - 2.*eta - 3))/256;
dN_deta{14} = -(81.*(xi - 1).*(3.*xi + 1).*(xi + 1).*(9.*eta.^2 - 2.*eta - 3))/256;
dN_deta{15} = (81.*(3.*xi + 1).*(xi - 1).*(xi + 1).*(9.*eta.^2 + 2.*eta - 3))/256;
dN_deta{16} = (81.*(9.*eta.^2 + 2.*eta - 3).*(- 3.*xi.^3 + xi.^2 + 3.*xi - 1))/256;

% Estas derivadas se calcularon con el codigo siguiente:
%{
syms xi eta

Ni1 = (1/16).*(xi-1).*(1-9.*xi.^2);
Ni2 = (9/16).*(1-xi.^2).*(1-3.*xi);
Ni3 = (9/16).*(1-xi.^2).*(1+3.*xi);
Ni4 = (1/16).*(xi+1).*(9.*xi.^2-1);

Nj1 = (1/16).*(eta-1).*(1-9.*eta.^2);
Nj2 = (9/16).*(1-eta.^2).*(1-3.*eta);
Nj3 = (9/16).*(1-eta.^2).*(1+3.*eta);
Nj4 = (1/16).*(eta+1).*(9.*eta.^2-1);

N = cell(16,1);
N{10} = Ni1.*Nj4;   N{9}  = Ni2.*Nj4;   N{8}  = Ni3.*Nj4;   N{7}  = Ni4.*Nj4;
N{11} = Ni1.*Nj3;   N{16} = Ni2.*Nj3;   N{15} = Ni3.*Nj3;   N{6}  = Ni4.*Nj3;
N{12} = Ni1.*Nj2;   N{13} = Ni2.*Nj2;   N{14} = Ni3.*Nj2;   N{5}  = Ni4.*Nj2;
N{1}  = Ni1.*Nj1;   N{2}  = Ni2.*Nj1;   N{3}  = Ni3.*Nj1;   N{4}  = Ni4.*Nj1;

dN_dxi  = cell(16,1);
dN_deta = cell(16,1);
for i=1:16
   dN_dxi{i}  = simple(diff(N{i},xi));
   fprintf('dN_dxi{%d} = %s;\n',i,char(dN_dxi{i}));
end

for i=1:16
   dN_deta{i} = simple(diff(N{i},eta));
   fprintf('dN_deta{%d} = %s;\n',i,char(dN_deta{i}));
end 
%}

dx_dxi  = zeros(n);
dy_dxi  = zeros(n);
dx_deta = zeros(n);
dy_deta = zeros(n);
for i = 1:16
   dx_dxi  = dx_dxi + dN_dxi{i}*xnod(i);
   dy_dxi  = dy_dxi + dN_dxi{i}*ynod(i);   
   dx_deta = dx_deta + dN_deta{i}*xnod(i);
   dy_deta = dy_deta + dN_deta{i}*ynod(i);   
end
% Calculo el determinante del Jacobiano
% J = [ dx_dxi   dy_dxi
%       dx_deta  dy_deta ]
detJ = [dx_dxi.*dy_deta - dx_deta.*dy_dxi];

% Se calcula el Jacobian ratio
JR = max(detJ(:))/min(detJ(:));
fprintf('JR = %f\n', JR);
if JR < 0 || JR > 40
    warning('Esta forma no es adecuada para un EF')
end

subplot(2,2,3);
pcolor(xi,eta,detJ)
colorbar
title(sprintf('Determinante de J. JR = %f', JR));
hold on
contour(xi,eta,detJ,[0 0], 'LineWidth',4, 'Color',[0 0 0]);
axis equal tight
