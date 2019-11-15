function N = N_H20(xi,eta,zeta)

N = [ ...
((xi - 1)*(xi - eta - eta*xi + eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 + eta*xi*zeta + 2))/8        % N1
-((xi^2 - 1)*(eta*zeta - eta + 1))/4                                                                 % N2
((xi + 1)*(eta + xi - eta*xi - eta*zeta^2 - eta^2*zeta + eta^2 + zeta^2 + eta*xi*zeta - 2))/8        % N3
((eta^2 - 1)*(xi + 1)*(zeta - 1))/4                                                                  % N4
((xi + 1)*(xi - eta + eta*xi + eta*zeta^2 - eta^2*zeta + eta^2 + zeta^2 - eta*xi*zeta - 2))/8        % N5
-((xi^2 - 1)*(eta - eta*zeta + 1))/4                                                                 % N6
((xi - 1)*(eta + xi + eta*xi - eta*zeta^2 + eta^2*zeta - eta^2 - zeta^2 - eta*xi*zeta + 2))/8        % N7
-((eta^2 - 1)*(xi - 1)*(zeta - 1))/4                                                                 % N8
-((zeta^2 - 1)*(eta - 1)*(xi - 1))/4                                                                 % N9
((zeta^2 - 1)*(eta - 1)*(xi + 1))/4                                                                  % N10
-((zeta^2 - 1)*(eta + 1)*(xi + 1))/4                                                                 % N11
((zeta^2 - 1)*(eta + 1)*(xi - 1))/4                                                                  % N12
-((xi - 1)*(eta - xi + eta*xi - eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 + eta*xi*zeta - 2))/8       % N13
((xi^2 - 1)*(eta + eta*zeta - 1))/4                                                                  % N14
((xi + 1)*(eta + xi - eta*xi - eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 - eta*xi*zeta - 2))/8        % N15
-((eta^2 - 1)*(xi + 1)*(zeta + 1))/4                                                                 % N16
((xi + 1)*(xi - eta + eta*xi + eta*zeta^2 + eta^2*zeta + eta^2 + zeta^2 + eta*xi*zeta - 2))/8        % N17
-((xi^2 - 1)*(eta + eta*zeta + 1))/4                                                                 % N18
((xi - 1)*(eta + xi + eta*xi - eta*zeta^2 - eta^2*zeta - eta^2 - zeta^2 + eta*xi*zeta + 2))/8        % N19
((eta^2 - 1)*(xi - 1)*(zeta + 1))/4 ];                                                               % N20

return;