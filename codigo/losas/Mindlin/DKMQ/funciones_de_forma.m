%% Funciones de forma principales del elemento finito
Nforma = @(xi,eta) [1/4*(1-xi)*(1-eta)        % N1
                    1/4*(1+xi)*(1-eta)        % N2
                    1/4*(1+xi)*(1+eta)        % N3
                    1/4*(1-xi)*(1+eta)];      % N4

%% Derivadas de N con respecto a xi    
dN_dxi = @(xi,eta) [-1/4*(1-eta)              % dN1_dxi
                     1/4*(1-eta)              % dN2_dxi
                     1/4*(1+eta)              % dN3_dxi
                    -1/4*(1+eta)    ];        % dN4_dxi
                        
%% Derivadas de N con respecto a eta    
dN_deta = @(xi,eta) [-1/4*(1-xi)              % dN1_deta
                     -1/4*(1+xi)              % dN2_deta
                      1/4*(1+xi)              % dN3_deta
                      1/4*(1-xi)    ];        % dN4_deta            

%% Funciones de forma secundarias del elemento finito
Pforma = @(xi,eta) [1/2*(1-xi^2)*(1-eta)      % P1
                    1/2*(1+xi)*(1-eta^2)      % P2
                    1/2*(1-xi^2)*(1+eta)      % P3
                    1/2*(1-xi)*(1-eta^2)];    % P4

%% Derivadas de P con respecto a xi
dP_dxi = @(xi,eta) [-xi*(1-eta)               % dP1_dxi
                    1/2*(1-eta^2)             % dP2_dxi
                    -xi*(1+eta)               % dP3_dxi
                    -1/2*(1-eta^2)  ];        % dP4_dxi
                        
%% Derivadas de P con respecto a eta    
dP_deta = @(xi,eta) [-1/2*(1-xi^2)            % dP1_deta
                     -eta*(1+xi)              % dP2_deta
                     1/2*(1-xi^2)             % dP3_deta
                     -eta*(1-xi)    ];        % dP4_deta       