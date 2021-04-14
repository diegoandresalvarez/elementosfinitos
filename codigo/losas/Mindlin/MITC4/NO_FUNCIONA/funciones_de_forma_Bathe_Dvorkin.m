%% Funciones de forma principales del elemento finito
Nforma = @(r,s) [  1/4*(1+r)*(1+s)       % = N3 Katili, et. al.
                   1/4*(1-r)*(1+s)       % = N4 Katili, et. al.
                   1/4*(1-r)*(1-s)       % = N1 Katili, et. al.
                   1/4*(1+r)*(1-s) ];    % = N2 Katili, et. al.
                    

%% Derivadas de N con respecto a r    
dN_dr = @(r,s) [  1/4*(1+s)              % dN1_dr
                 -1/4*(1+s)              % dN2_dr
                 -1/4*(1-s)              % dN3_dr
                  1/4*(1-s)    ];        % dN4_dr

                        
%% Derivadas de N con respecto a s    
dN_ds = @(r,s) [  1/4*(1+r)              % dN1_ds
                  1/4*(1-r)              % dN2_ds            
                 -1/4*(1-r)              % dN3_ds
                 -1/4*(1+r)    ];        % dN4_ds                      
