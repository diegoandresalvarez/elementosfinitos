% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa c5_funciones_forma_serendipito_rect_2D_8_nodos.m

%% NS son las funciones de forma del elemento serendípito de 8 nodos
NSforma = @(xi,eta) [ ...
-((eta - 1)*(xi - 1)*(eta + xi + 1))/4       % NS1
((xi^2 - 1)*(eta - 1))/2                     % NS2
((eta - 1)*(xi + 1)*(eta - xi + 1))/4        % NS3
-((eta^2 - 1)*(xi + 1))/2                    % NS4
((eta + 1)*(xi + 1)*(eta + xi - 1))/4        % NS5
-((xi^2 - 1)*(eta + 1))/2                    % NS6
((eta + 1)*(xi - 1)*(xi - eta + 1))/4        % NS7
((eta^2 - 1)*(xi - 1))/2                     % NS8
0                                        ];  % NS9 (inexistente)

%% Derivadas de NS con respecto a xi
dNS_dxi = @(xi,eta) [ ...
-((eta + 2*xi)*(eta - 1))/4                  % dNS1_dxi
eta*xi - xi                                  % dNS2_dxi
((eta - 2*xi)*(eta - 1))/4                   % dNS3_dxi
1/2 - eta^2/2                                % dNS4_dxi
((eta + 2*xi)*(eta + 1))/4                   % dNS5_dxi
-xi*(eta + 1)                                % dNS6_dxi
-((eta - 2*xi)*(eta + 1))/4                  % dNS7_dxi
eta^2/2 - 1/2                                % dNS8_dxi
0                                        ];  % dNS9_dxi (inexistente)

%% Derivadas de NS con respecto a eta
dNS_deta = @(xi,eta) [ ...
-((2*eta + xi)*(xi - 1))/4                   % dNS1_deta
xi^2/2 - 1/2                                 % dNS2_deta
((xi + 1)*(2*eta - xi))/4                    % dNS3_deta
-eta*(xi + 1)                                % dNS4_deta
((2*eta + xi)*(xi + 1))/4                    % dNS5_deta
1/2 - xi^2/2                                 % dNS6_deta
-((xi - 1)*(2*eta - xi))/4                   % dNS7_deta
eta*(xi - 1)                                 % dNS8_deta
0                                        ];  % dNS9_deta (inexistente)


% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa c9_funciones_forma_lagrangianos_9_nodos_rect_2D
%% NL son las funciones de forma del elemento lagrangiano de 9 nodos
NLforma = @(xi,eta) [ ...
(eta*xi*(eta - 1)*(xi - 1))/4               % NL1
-(eta*(xi^2 - 1)*(eta - 1))/2               % NL2
(eta*xi*(eta - 1)*(xi + 1))/4               % NL3
-(xi*(eta^2 - 1)*(xi + 1))/2                % NL4
(eta*xi*(eta + 1)*(xi + 1))/4               % NL5
-(eta*(xi^2 - 1)*(eta + 1))/2               % NL6
(eta*xi*(eta + 1)*(xi - 1))/4               % NL7
-(xi*(eta^2 - 1)*(xi - 1))/2                % NL8
(eta^2 - 1)*(xi^2 - 1)          ];          % NL9
 
%% Derivadas de NL con respecto a xi
dNL_dxi = @(xi,eta) [ ...
(eta*(2*xi - 1)*(eta - 1))/4                % dNL1_dxi
-eta*xi*(eta - 1)                           % dNL2_dxi
(eta*(2*xi + 1)*(eta - 1))/4                % dNL3_dxi
-((eta^2 - 1)*(2*xi + 1))/2                 % dNL4_dxi
(eta*(2*xi + 1)*(eta + 1))/4                % dNL5_dxi
-eta*xi*(eta + 1)                           % dNL6_dxi
(eta*(2*xi - 1)*(eta + 1))/4                % dNL7_dxi
-((eta^2 - 1)*(2*xi - 1))/2                 % dNL8_dxi
2*xi*(eta^2 - 1)                ];          % dNL9_dxi
 
%% Derivadas de NL con respecto a eta
dNL_deta = @(xi,eta) [ ...
(xi*(2*eta - 1)*(xi - 1))/4                 % dNL1_deta
-((2*eta - 1)*(xi^2 - 1))/2                 % dNL2_deta
(xi*(2*eta - 1)*(xi + 1))/4                 % dNL3_deta
-eta*xi*(xi + 1)                            % dNL4_deta
(xi*(2*eta + 1)*(xi + 1))/4                 % dNL5_deta
-((2*eta + 1)*(xi^2 - 1))/2                 % dNL6_deta
(xi*(2*eta + 1)*(xi - 1))/4                 % dNL7_deta
-eta*xi*(xi - 1)                            % dNL8_deta
2*eta*(xi^2 - 1)                ];          % dNL9_deta
