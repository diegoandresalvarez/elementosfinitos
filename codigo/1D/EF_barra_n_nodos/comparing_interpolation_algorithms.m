% Interpolation algorithms already implemented in MATLAB
%
% WHO   DATE            WHAT
% DAA   Mar 1,  2010    First algorithm
% DAA   Nov 11, 2021    Introducing MAKIMA
%
% DAA - Diego Andres Alvarez Marin - daalvarez@unal.edu.co

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
end

xi = -10:0.1:10;
yi = interp1(x,y,xi,'nearest','extrap'); h1 = plot(xi,yi,'b.', 'LineWidth',2);
yi = interp1(x,y,xi,'linear', 'extrap'); h2 = plot(xi,yi,'r-', 'LineWidth',2);
yi = interp1(x,y,xi,'spline', 'extrap'); h3 = plot(xi,yi,'m-', 'LineWidth',2);
yi = interp1(x,y,xi,'makima', 'extrap'); h4 = plot(xi,yi,'r--','LineWidth',2);
yi = interp1(x,y,xi,'pchip',  'extrap'); h5 = plot(xi,yi,'g-', 'LineWidth',2);
% The 'cubic' method does not support extrapolation.
yi = interp1(x,y,xi,'cubic'           ); h6 = plot(xi,yi,'c-', 'LineWidth',2);

% Polynomial curve fitting (LAGRANGE POLYNOMIAL)
% Calculate the coefficients in the approximating polynomial of degree n-1
p = polyfit(x,y,n-1);
h7 = plot(xi,polyval(p,xi),'k-','LineWidth',4);
axis([-10  10  -5  5]);

legend([h1 h2 h3 h4 h5 h6 h7], ...
    'nearest','linear','spline','makima','pchip','cubic','Lagrange', ...
    'Location','NorthEastOutside')

% TAREA: leer la ayuda del los comandos pchip, makima, spline, interp1
