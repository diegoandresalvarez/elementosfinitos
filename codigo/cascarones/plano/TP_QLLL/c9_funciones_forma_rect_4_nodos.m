%% N son las funciones de forma del elemento rectangular de 4 nodos
Nforma = @(xi,eta) [  ((xi - 1)*(eta - 1))/4    % N1
                     -((xi + 1)*(eta - 1))/4    % N2
                      ((xi + 1)*(eta + 1))/4    % N3
                     -((xi - 1)*(eta + 1))/4 ]; % N4

%% Derivadas de N con respecto a xi
dN_dxi  = @(xi,eta) [ (eta - 1)/4               % dN1_dxi
                     -(eta - 1)/4               % dN2_dxi
                      (eta + 1)/4               % dN3_dxi
                     -(eta + 1)/4 ];            % dN4_dxi
 
%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [  (xi - 1)/4               % dN1_deta
                      -(xi + 1)/4               % dN2_deta
                       (xi + 1)/4               % dN3_deta
                      -(xi - 1)/4 ];            % dN4_deta
                      
%{
syms xi eta
L2_xi_1  = factor(poly2sym(polyfit([-1 1],[1 0],1),xi));
L2_xi_2  = factor(poly2sym(polyfit([-1 1],[0 1],1),xi));
L2_eta_1 = factor(poly2sym(polyfit([-1 1],[1 0],1),eta));
L2_eta_2 = factor(poly2sym(polyfit([-1 1],[0 1],1),eta));

N = sym(zeros(1,4));
N(1) = L2_xi_1 * L2_eta_1;
N(2) = L2_xi_2 * L2_eta_1;
N(3) = L2_xi_2 * L2_eta_2;
N(4) = L2_xi_1 * L2_eta_2;

dN_dxi  = diff(N,xi);
dN_deta = diff(N,eta);
%}
