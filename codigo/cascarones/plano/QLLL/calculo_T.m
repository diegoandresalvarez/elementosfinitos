function [T, lambdae] = calculo_T(xnod_e)
%% Calcula la matriz de transformacion de coordenadas (sec 10.7.1)
%
% [T, lambdae] = calculo_T(xnod(LaG(e,:),:);

%%
vij = xnod_e(2,:) - xnod_e(1,:);
vim = xnod_e(3,:) - xnod_e(1,:);
vxp = vij/norm(vij);   
vzp = cross(vij, vim); vzp = vzp/norm(vzp);
vyp = cross(vzp,vxp);

% el sparse() es para obligar a que Le y T sean sparse
lambdae  = sparse([  vxp; vyp; vzp ]);
lambdage = sparse([ -vyp; vxp ]);
   
Le = blkdiag(lambdae, lambdage); 
T  = blkdiag(Le,Le,Le,Le);               % 4 veces para los 4 nodos
return