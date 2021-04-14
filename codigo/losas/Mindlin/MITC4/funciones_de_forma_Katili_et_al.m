%% Funciones de forma principales del elemento finito
Nforma = @(xi,eta) [  1/4*(1-xi)*(1-eta)        % N1
                      1/4*(1+xi)*(1-eta)        % N2
                      1/4*(1+xi)*(1+eta)        % N3
                      1/4*(1-xi)*(1+eta) ];     % N4

%% Derivadas de N con respecto a xi    
dN_dxi = @(xi,eta) [ -1/4*(1-eta)              % dN1_dxi
                      1/4*(1-eta)              % dN2_dxi
                      1/4*(1+eta)              % dN3_dxi
                     -1/4*(1+eta)    ];        % dN4_dxi
                        
%% Derivadas de N con respecto a eta    
dN_deta = @(xi,eta) [ -1/4*(1-xi)              % dN1_deta
                      -1/4*(1+xi)              % dN2_deta
                       1/4*(1+xi)              % dN3_deta
                       1/4*(1-xi)    ];        % dN4_deta            
