function [T, lambdae] = calculo_T(xnod_e)
%% Calcula la matriz de transformación de coordenadas (sec 10.7.1)
%
% [T, lambdae] = calculo_T(xnod(LaG(e,:),:);

%%
vij = xnod_e(2,:) - xnod_e(1,:);
vim = xnod_e(3,:) - xnod_e(1,:);
vxp = vij/norm(vij);   
vzp = cross(vij, vim); vzp = vzp/norm(vzp);
vyp = cross(vzp,vxp);
lambdae  = [ vxp; vyp; vzp ];
lambdage = [ -vyp; vxp ];
   
Le = blkdiag(lambdae, sparse(lambdage)); % forzamos a que T sea sparse
T  = blkdiag(Le,Le,Le,Le);               % 4 veces para los 4 nodos

return