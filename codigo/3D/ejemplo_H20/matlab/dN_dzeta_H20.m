function dN_dzeta = dN_dzeta_H20(xi,eta,zeta)
% Derivadas de las funciones de forma con respecto a la variable "zeta" del
% EF hexaedrico serendipito isoparametrico de 20 nodos

dN_dzeta = [ ...
  ((eta - 1)*(xi - 1)*(eta + xi + 2*zeta + 1))/8    % dN1_dzeta 
 -((xi^2 - 1)*(eta - 1))/4                          % dN2_dzeta 
 -((eta - 1)*(xi + 1)*(eta - xi + 2*zeta + 1))/8    % dN3_dzeta 
  ((eta^2 - 1)*(xi + 1))/4                          % dN4_dzeta 
 -((eta + 1)*(xi + 1)*(eta + xi - 2*zeta - 1))/8    % dN5_dzeta 
  ((xi^2 - 1)*(eta + 1))/4                          % dN6_dzeta 
 -((eta + 1)*(xi - 1)*(xi - eta + 2*zeta + 1))/8    % dN7_dzeta 
 -((eta^2 - 1)*(xi - 1))/4                          % dN8_dzeta 
 -(zeta*(eta - 1)*(xi - 1))/2                       % dN9_dzeta 
  (zeta*(eta - 1)*(xi + 1))/2                       % dN10_dzeta
 -(zeta*(eta + 1)*(xi + 1))/2                       % dN11_dzeta
  (zeta*(eta + 1)*(xi - 1))/2                       % dN12_dzeta
 -((eta - 1)*(xi - 1)*(eta + xi - 2*zeta + 1))/8    % dN13_dzeta
  ((xi^2 - 1)*(eta - 1))/4                          % dN14_dzeta
  ((eta - 1)*(xi + 1)*(eta - xi - 2*zeta + 1))/8    % dN15_dzeta
 -((eta^2 - 1)*(xi + 1))/4                          % dN16_dzeta
  ((eta + 1)*(xi + 1)*(eta + xi + 2*zeta - 1))/8    % dN17_dzeta
 -((xi^2 - 1)*(eta + 1))/4                          % dN18_dzeta
 -((eta + 1)*(xi - 1)*(eta - xi + 2*zeta - 1))/8    % dN19_dzeta
  ((eta^2 - 1)*(xi - 1))/4                          % dN20_dzeta
];

return;
