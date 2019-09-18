function [xi,w,P] = gausslegendre_quad(m)
% Integration using a Gauss-Legendre quadrature:
%
% Run the following example:
%{
m = 4;

[xi,w,P] = gausslegendre_quad(m);

x = linspace(-1,1,100);
y = zeros(m+1,100);
for i = 1:m+1
   y(i,:) = polyval(P{i},x);
end;

figure
h = plot(x, y);
hold on;
plot(xi,zeros(size(xi)),'ro');
legend(h, num2str((0:m)'), 'Location', 'Best');
title(sprintf('Primeros %d polinomios de Legendre', m+1))
grid on

axis([-1 1 -1.1 1.1]);

a = 0; b = 0.8;
f = @(x) 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5;
sol = 3076/1875;
fprintf('Error = %g\n', abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - sol))

a = 0; b = pi/2;
f = @(x) sin(x);
sol = 1;
fprintf('Error = %g\n', abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - sol))
%}

% WHO   DATE            WHAT
% DAAM   Mar 11, 2010    First algorithm
% DAAM   Sep 18, 2019    Better readability, MATLAB 2019a
%
% DAAM - Diego Andres Alvarez Marin - daalvarez@unal.edu.co

%% Gaussian quadrature xi and w
%% Polynomials (Bonnet's recursion): 
P = cell(m+1, 1);
P{0 + 1} = 1;
P{1 + 1} = [1 0];
for n = 2:m
   P{n + 1} = ((2*n-1)*[P{n-1 + 1} 0] - (n-1)*[0 0 P{n-2 + 1}])/n;
end

%% Roots
xi = sort(roots(P{m + 1}));

%% Weights: VERSION 1
s = polyder(P{m + 1});
w = 2./((1 - xi.^2).*polyval(s,xi).^2);

%% Weights: VERSION 2
%{
A = zeros(m,m);
b = zeros(m,1);
for i=1:m
  A(i,:) = xi.^(i-1)';
  b(i) = (1-(-1)^i)/i;
end
w = A\b;
%}

if ~isreal(w)
   error('m is too large. The weights cannot be complex')
end

return;

%% Version que utiliza el toolbox de algebra simbolica del MATLAB
%{
m = 5
if m >= 6
   error('m is too large. The symbolic toolbox does not work with this configuration')
end

syms xi;

% Polynomials (Bonnet's recursion): 
P = cell(m+1,1);
P{0 + 1} = 1;
P{1 + 1} = xi;
for n = 2:m
   P{n + 1} = (1/n)*((2*n-1)*xi*P{n-1 + 1} - (n-1)*P{n-2 + 1});
end

% Roots
xir = simplify(solve(P{m + 1} == 0, xi))

% Weights
w = simplify(2./((1-xir.^2).*subs(diff(P{m + 1}), 'xi', xir).^2));

% Print the results
[xir w]
%}
