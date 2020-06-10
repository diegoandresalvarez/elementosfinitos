% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa "codigo/2D/deduccion_funciones_forma/FF_lagrangianos_Q9.m"

%% N son las funciones de forma del elemento lagrangiano de 9 nodos
Nforma = @(xi,eta) [ ...
        (eta*xi*(eta - 1)*(xi - 1))/4               % N1
        -(eta*(xi^2 - 1)*(eta - 1))/2               % N2
        (eta*xi*(eta - 1)*(xi + 1))/4               % N3
        -(xi*(eta^2 - 1)*(xi + 1))/2                % N4
        (eta*xi*(eta + 1)*(xi + 1))/4               % N5
        -(eta*(xi^2 - 1)*(eta + 1))/2               % N6
        (eta*xi*(eta + 1)*(xi - 1))/4               % N7
        -(xi*(eta^2 - 1)*(xi - 1))/2                % N8
        (eta^2 - 1)*(xi^2 - 1)          ];          % N9
 
%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
        (eta*(2*xi - 1)*(eta - 1))/4                % dN1_dxi
        -eta*xi*(eta - 1)                           % dN2_dxi
        (eta*(2*xi + 1)*(eta - 1))/4                % dN3_dxi
        -((eta^2 - 1)*(2*xi + 1))/2                 % dN4_dxi
        (eta*(2*xi + 1)*(eta + 1))/4                % dN5_dxi
        -eta*xi*(eta + 1)                           % dN6_dxi
        (eta*(2*xi - 1)*(eta + 1))/4                % dN7_dxi
        -((eta^2 - 1)*(2*xi - 1))/2                 % dN8_dxi
        2*xi*(eta^2 - 1)                ];          % dN9_dxi
 
%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
        (xi*(2*eta - 1)*(xi - 1))/4                 % dN1_deta
        -((2*eta - 1)*(xi^2 - 1))/2                 % dN2_deta
        (xi*(2*eta - 1)*(xi + 1))/4                 % dN3_deta
        -eta*xi*(xi + 1)                            % dN4_deta
        (xi*(2*eta + 1)*(xi + 1))/4                 % dN5_deta
        -((2*eta + 1)*(xi^2 - 1))/2                 % dN6_deta
        (xi*(2*eta + 1)*(xi - 1))/4                 % dN7_deta
        -eta*xi*(xi - 1)                            % dN8_deta
        2*eta*(xi^2 - 1)                ];          % dN9_deta
