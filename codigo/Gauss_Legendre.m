function [z,w] = Gauss_Legendre(N)
%% Calcula las raices del polinomio y los pesos de la cuadratura de Gauss Legendre
%
%[z,w] = Gauss_Legendre(N)
%
%z		Las raices del polinomio
%w		Los pesos de Gauss - Legendre

if N == 1
   z = 0;
   w = 2;
   return;
end;

[c,s] = legendre(N)							%El polinomio de Legendre
z = sort(roots(c));							%Las raices del polinomio
w = 2./((1-z.^2).*polyval(s,z).^2);		%Los pesos de Gauss - Legendre

return;

%% -----------------------------------------------------------------------------

function [PN,DN] = legendre(N)
%% Calcula los coeficientes del polinomio de Legendre de orden N y su derivada
%
%[PN,DN]=legendre(N)
%
%N	Orden del polinomio
%PN	Polinomio
%DN	Derivada de PN

if N == 1
   PN = 1;
   DN = 0;
   return;
end;

Pnm1 = 1;
Pn = [1 0];

for n = 1:N-1
   Pnp1 = ((2*n+1)*[Pn 0]-n*[0 0 Pnm1])/(n+1);
   Pnm1 = Pn;
   Pn = Pnp1;
end

PN = Pnp1;

DN = zeros(1,N);
for n = 1:N
   DN(n) = PN(n)*(N-n+1);
end

return;

