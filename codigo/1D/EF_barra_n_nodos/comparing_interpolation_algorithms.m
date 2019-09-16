% Interpolation algorithms already implemented in MATLAB
%
% WHO   DATE            WHAT
% DAA   Mar 1, 2010    First algorithm
%
% DAA - Diego Andres Alvarez Marin - diegoandresalvarez@gmx.net

% Input the data points
n = 10; %number of data points
x = zeros(n,1);
y = zeros(n,1);
figure
title(sprintf('Click with the mouse %d times. The interpolation will pass through those points',n))
axis([-5 5 -5 5]);
hold on;
for i = 1:n
   [x(i),y(i)] = ginput(1);
   plot(x(i), y(i), 'ko','MarkerSize',12);  
   plot(x(i), y(i), 'kx','MarkerSize',12);
end;

xi = (-10: 0.1: 10)';
yi = interp1(x,y,xi,'nearest','extrap'); h1  = plot(xi,yi,'b.','LineWidth',2);
yi = interp1(x,y,xi,'linear','extrap');  h2  = plot(xi,yi,'r-','LineWidth',2);
yi = interp1(x,y,xi,'spline','extrap');  h3a = plot(xi,yi,'m-','LineWidth',2);
yi = spline(x,y,xi);                     h3b = plot(xi,yi,'r--','LineWidth',2);
yi = interp1(x,y,xi,'pchip','extrap');   h4  = plot(xi,yi,'g-','LineWidth',2);
yi = interp1(x,y,xi,'cubic','extrap');   h5  = plot(xi,yi,'c-','LineWidth',2);

% Polynomial curve fitting (LAGRANGE POLYNOMIAL)
% Calculate the coefficients in the approximating polynomial of degree n-1
p = polyfit(x,y,n-1);
h6 = plot(xi,polyval(p,xi),'k-','LineWidth',4);
axis([-10  10  -5  5]);

legend([h1 h2 h3a h3b h4 h5 h6], ...
    'nearest','linear','spline','spline function','pchip','cubic','Lagrange', ...
    'Location','NorthEastOutside')

% TAREA: leer la ayuda del los comandos pchip, spline, interp1