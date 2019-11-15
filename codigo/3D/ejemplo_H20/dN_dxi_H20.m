function dN_dxi = dN_dxi_H20(xi,eta,zeta)

% Derivadas con respecto a xi:
dN_dxi = [ ...
(2*xi - 2*eta*xi - eta*zeta + eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 + 2*eta*xi*zeta + 1)/8        % dN1_dxi
-(xi*(eta*zeta - eta + 1))/2                                                                         % dN2_dxi
-(2*eta*xi - 2*xi - eta*zeta + eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 - 2*eta*xi*zeta + 1)/8       % dN3_dxi
((eta^2 - 1)*(zeta - 1))/4                                                                           % dN4_dxi
-(eta*zeta - 2*eta*xi - 2*xi - eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 + 2*eta*xi*zeta + 1)/8       % dN5_dxi
(xi*(eta*zeta - eta - 1))/2                                                                          % dN6_dxi
(2*xi + 2*eta*xi + eta*zeta - eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 - 2*eta*xi*zeta + 1)/8        % dN7_dxi
-((eta^2 - 1)*(zeta - 1))/4                                                                          % dN8_dxi
-((zeta^2 - 1)*(eta - 1))/4                                                                          % dN9_dxi
((zeta^2 - 1)*(eta - 1))/4                                                                           % dN10_dxi
-((zeta^2 - 1)*(eta + 1))/4                                                                          % dN11_dxi
((zeta^2 - 1)*(eta + 1))/4                                                                           % dN12_dxi
-(2*eta*xi - 2*xi - eta*zeta - eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 + 2*eta*xi*zeta - 1)/8       % dN13_dxi
(xi*(eta + eta*zeta - 1))/2                                                                          % dN14_dxi
(2*xi - 2*eta*xi - eta*zeta - eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 - 2*eta*xi*zeta - 1)/8        % dN15_dxi
-((eta^2 - 1)*(zeta + 1))/4                                                                          % dN16_dxi
(2*xi + 2*eta*xi + eta*zeta + eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 + 2*eta*xi*zeta - 1)/8        % dN17_dxi
-(xi*(eta + eta*zeta + 1))/2                                                                         % dN18_dxi
-(eta*zeta - 2*eta*xi - 2*xi + eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 - 2*eta*xi*zeta - 1)/8       % dN19_dxi
((eta^2 - 1)*(zeta + 1))/4 ];                                                                        % dN20_dxi

return;