% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa c5_funciones_forma_serendipito_rect_2D_8_nodos.m
global Nforma dN_dxi dN_deta

%% N son las funciones de forma del elemento serendipito de 8 nodos
Nforma = @(xi,eta) [ ...
-((eta - 1)*(xi - 1)*(eta + xi + 1))/4           % N1
((xi^2 - 1)*(eta - 1))/2                         % N2
((eta - 1)*(xi + 1)*(eta - xi + 1))/4            % N3
-((eta^2 - 1)*(xi + 1))/2                        % N4
((eta + 1)*(xi + 1)*(eta + xi - 1))/4            % N5
-((xi^2 - 1)*(eta + 1))/2                        % N6
((eta + 1)*(xi - 1)*(xi - eta + 1))/4            % N7
((eta^2 - 1)*(xi - 1))/2                ];       % N8
 
%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
-((eta + 2*xi)*(eta - 1))/4                      % dN1_dxi
eta*xi - xi                                      % dN2_dxi
((eta - 2*xi)*(eta - 1))/4                       % dN3_dxi
1/2 - eta^2/2                                    % dN4_dxi
((eta + 2*xi)*(eta + 1))/4                       % dN5_dxi
-xi*(eta + 1)                                    % dN6_dxi
-((eta - 2*xi)*(eta + 1))/4                      % dN7_dxi
eta^2/2 - 1/2                           ];       % dN8_dxi

%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
 -((2*eta + xi)*(xi - 1))/4                      % dN1_deta
 xi^2/2 - 1/2                                    % dN2_deta
 ((xi + 1)*(2*eta - xi))/4                       % dN3_deta
 -eta*(xi + 1)                                   % dN4_deta
 ((2*eta + xi)*(xi + 1))/4                       % dN5_deta
 1/2 - xi^2/2                                    % dN6_deta
 -((xi - 1)*(2*eta - xi))/4                      % dN7_deta
 eta*(xi - 1)                           ];       % dN8_deta
