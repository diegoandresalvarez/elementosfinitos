function [xi,w,P] = gausslegendre_quad(m)
% Integration using a Gauss-Legendre quadrature:
%
% Run the following example:
%{
m = 4;

[xi,w,P] = gausslegendre_quad(m);

xx = linspace(-1,1,100);
yy = zeros(m+1,100);
for i = 1:m+1
   yy(i,:) = polyval(P{i},xx);
end;

figure
plot(xx,yy);
hold on;
legend(num2str((0:m)'));
plot(xi,zeros(size(xi)),'ro');
grid on

axis([-1 1 -1.1 1.1]);

a = 0; b = 0.8;
f = @(x) 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5;
abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - 3076/1875)

a = 0; b = pi/2;
f = @(x) sin(x);
abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - 1)
%}

% WHO   DATE            WHAT
% DAA   Mar 11, 2010    First algorithm
%
% DAA - Diego Andres Alvarez Marin - daalvarez@unal.edu.co

%% Gaussian quadrature xi and w
%% Polynomials (Bonnet's recursion): 
% Observe that P_0 = P{1}, P_1 = P{2}, ... P_{n-1} = P{n}, P_n = P{n+1}
% in as much as MATLAB does not make 0-based indexing of arrays:
P = cell(m+1,1);
P{1} = 1;
P{2} = [1 0];
for n = 2:m
   P{n+1} = ((2*n-1)*[P{n} 0] - (n-1)*[0 0 P{n-1}])/n;
end

%% Roots
xi = sort(roots(P{m+1}));

%% Weights: VERSION 1
s = polyder(P{m+1});
w = 2./((1-xi.^2).*polyval(s,xi).^2);

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

P = cell(m+1,1);
P{1} = 1;
P{2} = xi;
for n = 2:m
   P{n+1} = (1/n)*((2*n-1)*xi*P{n} - (n-1)*P{n-1});
end

xir = simple(solve(P{m+1},'xi'));

w = simple(2./((1-xir.^2).*subs(diff(P{m+1}), 'xi', xir).^2));

[xir w]
%}
%%

