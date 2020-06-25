function [T, lambdae]  = calculo_T2(xnod_e)
%% Calcula la matriz de transformación de coordenadas (sec 10.7.2)
%
% [T, lambdae] = calculo_T2(xnod(LaG(e,:),:);

%%
v12 = xnod_e(2,:) - xnod_e(1,:);
v13 = xnod_e(3,:) - xnod_e(1,:);

vzp = cross(v12,v13); vzp = vzp/norm(vzp);

vxp = [ 1/sqrt(1+(vzp(1)/vzp(3))^2), 0, -1/sqrt(1+(vzp(3)/vzp(1))^2) ];

dd = vxp(1)*vzp(1) + vxp(3)*vzp(3);
if abs(dd) > 1e-8
   vxp(3) = -vxp(3);
end

if (vzp(3) == 0) && (vzp(1) == 0)
   vxp = [1 0 0];
end

vyp = cross(vzp,vxp); vyp = vyp/norm(vyp);

if vyp(2) < 0
   vyp = -vyp;
   vxp = -vxp;
end

lambdae  = [ vxp; vyp; vzp ];
lambdage = [ -vyp; vxp ];
   
Le = blkdiag(lambdae, sparse(lambdage)); % forzamos a que T sea sparse
T  = blkdiag(Le,Le,Le,Le);               % 4 veces para los 4 nodos

return