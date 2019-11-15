function dN_dzeta = dN_dzeta_H20(xi,eta,zeta)

% Derivadas con respecto a zeta:
dN_dzeta = [ ...
((xi - 1)*(eta*xi - 2*zeta + 2*eta*zeta + eta^2))/8   % dN1_dzeta
-(eta*(xi^2 - 1))/4                                   % dN2_dzeta
((xi + 1)*(2*zeta + eta*xi - 2*eta*zeta - eta^2))/8   % dN3_dzeta
((eta^2 - 1)*(xi + 1))/4                              % dN4_dzeta
((xi + 1)*(2*zeta - eta*xi + 2*eta*zeta - eta^2))/8   % dN5_dzeta
(eta*(xi^2 - 1))/4                                    % dN6_dzeta
((xi - 1)*(eta^2 - eta*xi - 2*eta*zeta - 2*zeta))/8   % dN7_dzeta
-((eta^2 - 1)*(xi - 1))/4                             % dN8_dzeta
-(zeta*(eta - 1)*(xi - 1))/2                          % dN9_dzeta
(zeta*(eta - 1)*(xi + 1))/2                           % dN10_dzeta
-(zeta*(eta + 1)*(xi + 1))/2                          % dN11_dzeta
(zeta*(eta + 1)*(xi - 1))/2                           % dN12_dzeta
-((xi - 1)*(2*zeta + eta*xi - 2*eta*zeta + eta^2))/8  % dN13_dzeta
(eta*(xi^2 - 1))/4                                    % dN14_dzeta
((xi + 1)*(2*zeta - eta*xi - 2*eta*zeta + eta^2))/8   % dN15_dzeta
-((eta^2 - 1)*(xi + 1))/4                             % dN16_dzeta
((xi + 1)*(2*zeta + eta*xi + 2*eta*zeta + eta^2))/8   % dN17_dzeta
-(eta*(xi^2 - 1))/4                                   % dN18_dzeta
-((xi - 1)*(2*zeta - eta*xi + 2*eta*zeta + eta^2))/8  % dN19_dzeta
((eta^2 - 1)*(xi - 1))/4 ];                           % dN20_dzeta

return;