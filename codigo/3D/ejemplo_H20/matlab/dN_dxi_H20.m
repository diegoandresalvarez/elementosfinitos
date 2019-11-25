function dN_dxi = dN_dxi_H20(xi,eta,zeta)
% Derivadas de las funciones de forma con respecto a la variable "xi" del
% EF hexaedrico serendipito isoparametrico de 20 nodos

dN_dxi = [ ...
  ((eta - 1)*(zeta - 1)*(eta + 2*xi + zeta + 1))/8   % dN1_dxi 
 -(xi*(eta - 1)*(zeta - 1))/2                        % dN2_dxi 
 -((eta - 1)*(zeta - 1)*(eta - 2*xi + zeta + 1))/8   % dN3_dxi 
  ((eta^2 - 1)*(zeta - 1))/4                         % dN4_dxi 
 -((eta + 1)*(zeta - 1)*(eta + 2*xi - zeta - 1))/8   % dN5_dxi 
  (xi*(eta + 1)*(zeta - 1))/2                        % dN6_dxi 
 -((eta + 1)*(zeta - 1)*(2*xi - eta + zeta + 1))/8   % dN7_dxi 
 -((eta^2 - 1)*(zeta - 1))/4                         % dN8_dxi 
 -((zeta^2 - 1)*(eta - 1))/4                         % dN9_dxi 
  ((zeta^2 - 1)*(eta - 1))/4                         % dN10_dxi
 -((zeta^2 - 1)*(eta + 1))/4                         % dN11_dxi
  ((zeta^2 - 1)*(eta + 1))/4                         % dN12_dxi
 -((eta - 1)*(zeta + 1)*(eta + 2*xi - zeta + 1))/8   % dN13_dxi
  (xi*(eta - 1)*(zeta + 1))/2                        % dN14_dxi
  ((eta - 1)*(zeta + 1)*(eta - 2*xi - zeta + 1))/8   % dN15_dxi
 -((eta^2 - 1)*(zeta + 1))/4                         % dN16_dxi
  ((eta + 1)*(zeta + 1)*(eta + 2*xi + zeta - 1))/8   % dN17_dxi
 -(xi*(eta + 1)*(zeta + 1))/2                        % dN18_dxi
 -((eta + 1)*(zeta + 1)*(eta - 2*xi + zeta - 1))/8   % dN19_dxi
  ((eta^2 - 1)*(zeta + 1))/4                         % dN20_dxi
];

return;
