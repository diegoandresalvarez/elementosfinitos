function dN_deta = dN_deta_H20(xi,eta,zeta)

% Derivadas con respecto a eta:
dN_deta = [ ...
((xi - 1)*(zeta - 1)*(2*eta + xi + zeta + 1))/8       % dN1_deta
-((xi^2 - 1)*(zeta - 1))/4                            % dN2_deta
-((xi + 1)*(zeta - 1)*(2*eta - xi + zeta + 1))/8      % dN3_deta
(eta*(xi + 1)*(zeta - 1))/2                           % dN4_deta
-((xi + 1)*(zeta - 1)*(2*eta + xi - zeta - 1))/8      % dN5_deta
((xi^2 - 1)*(zeta - 1))/4                             % dN6_deta
((xi - 1)*(zeta - 1)*(2*eta - xi - zeta - 1))/8       % dN7_deta
-(eta*(xi - 1)*(zeta - 1))/2                          % dN8_deta
-((zeta^2 - 1)*(xi - 1))/4                            % dN9_deta
((zeta^2 - 1)*(xi + 1))/4                             % dN10_deta
-((zeta^2 - 1)*(xi + 1))/4                            % dN11_deta
((zeta^2 - 1)*(xi - 1))/4                             % dN12_deta
-((xi - 1)*(zeta + 1)*(2*eta + xi - zeta + 1))/8      % dN13_deta
((xi^2 - 1)*(zeta + 1))/4                             % dN14_deta
((xi + 1)*(zeta + 1)*(2*eta - xi - zeta + 1))/8       % dN15_deta
-(eta*(xi + 1)*(zeta + 1))/2                          % dN16_deta
((xi + 1)*(zeta + 1)*(2*eta + xi + zeta - 1))/8       % dN17_deta
-((xi^2 - 1)*(zeta + 1))/4                            % dN18_deta
-((xi - 1)*(zeta + 1)*(2*eta - xi + zeta - 1))/8      % dN19_deta
(eta*(xi - 1)*(zeta + 1))/2 ];                        % dN20_deta

return;