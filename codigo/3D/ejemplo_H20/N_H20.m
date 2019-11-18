function N = N_H20(xi,eta,zeta)
% Funciones de forma del EF hexaedrico serendipito isoparametrico de 20 nodos

N = [ ...
     ((eta - 1)*(xi - 1)*(zeta - 1)*(eta + xi + zeta + 2))/8     % N1
    -((xi^2 - 1)*(eta - 1)*(zeta - 1))/4                         % N2
    -((eta - 1)*(xi + 1)*(zeta - 1)*(eta - xi + zeta + 2))/8     % N3
     ((eta^2 - 1)*(xi + 1)*(zeta - 1))/4                         % N4
    -((eta + 1)*(xi + 1)*(zeta - 1)*(eta + xi - zeta - 2))/8     % N5
     ((xi^2 - 1)*(eta + 1)*(zeta - 1))/4                         % N6
    -((eta + 1)*(xi - 1)*(zeta - 1)*(xi - eta + zeta + 2))/8     % N7
    -((eta^2 - 1)*(xi - 1)*(zeta - 1))/4                         % N8
    -((zeta^2 - 1)*(eta - 1)*(xi - 1))/4                         % N9
     ((zeta^2 - 1)*(eta - 1)*(xi + 1))/4                         % N10
    -((zeta^2 - 1)*(eta + 1)*(xi + 1))/4                         % N11
     ((zeta^2 - 1)*(eta + 1)*(xi - 1))/4                         % N12
    -((eta - 1)*(xi - 1)*(zeta + 1)*(eta + xi - zeta + 2))/8     % N13
     ((xi^2 - 1)*(eta - 1)*(zeta + 1))/4                         % N14
     ((eta - 1)*(xi + 1)*(zeta + 1)*(eta - xi - zeta + 2))/8     % N15
    -((eta^2 - 1)*(xi + 1)*(zeta + 1))/4                         % N16
     ((eta + 1)*(xi + 1)*(zeta + 1)*(eta + xi + zeta - 2))/8     % N17
    -((xi^2 - 1)*(eta + 1)*(zeta + 1))/4                         % N18
    -((eta + 1)*(xi - 1)*(zeta + 1)*(eta - xi + zeta - 2))/8     % N19
     ((eta^2 - 1)*(xi - 1)*(zeta + 1))/4                         % N20
];                                                              

return;
