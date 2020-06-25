function [vg1, vg2, vg3, t] = vg123_t(rlo, rup)
%% Calcula los vectores vg1, vg2 y vg3 para cada nodo, las coord de los
%% nodos del plano de referencia y el espesor del cascaron
%
% [vg1, vg2, vg3, t] = vg123_t(rlo, rup)
% 
% rlo, rup  coordenadas de los nodos que definen las caras de los EFs
% vgi[nno x 3] i = 1,2,3 base que definen los ejes del sist de coord nodales
% t         espesor del cascaron en cada uno de los nodos

%%
try
   V3   = rup - rlo;
catch err
   error('rup y rlo deben tener el mismo tamano')
end
   
nno = size(rup, 1);

vg1 = zeros(nno, 3);
vg2 = zeros(nno, 3);
vg3 = zeros(nno, 3);
t   = zeros(nno, 1);
for i = 1:nno
   t(i) = norm(V3(i,:));               % se calcula el espesor en el nodo i
   vg3(i,:) = V3(i,:)/t(i);            % se normaliza vg3
   
   % si vg3(i,:) == [0 1 0] o vg3(i,:) == [0 -1 0]
   if (norm(vg3(i,:) - [0 1 0]) < 1e-6) || (norm(vg3(i,:) - [0 -1 0]) < 1e-6)
      vg1(i,:) = [1 0 0];
   else
      tmp = cross([0 1 0], vg3(i,:));
      vg1(i,:) = tmp/norm(tmp);
   end
   vg2(i,:) = cross(vg3(i,:), vg1(i,:));
end
      
return
