function [xnod,LaG] = import_from_GiD(filename)
% Esta rutina importa los archivos de malla generados por GiD
%
% Uso:
% [xnod,LaG] = import_from_GiD(filename)
%
% POR EL MOMENTO SOLO FUNCIONA PARA ELEMENTOS H20

% VER EN LA AYUDA DE GID:
% Postprocess mesh format: ProjectName.post.msh
% Postprocess results format: ProjectName.post.res
% http://www-opale.inrialpes.fr/Aerochina/info/en/html-version/gid_17.html

%%
fid = fopen([filename '.msh'],'r');

linea = fgetl(fid);  % MESH dimension 3 ElemType Hexahedra Nnode 20 Coordinates
expr = 'MESH dimension (\d+) ElemType (\w+) Nnode (\d+)';

tmp = regexp(linea, expr, 'tokens');
[ndim, ElemType, nno_por_elem] = tmp{1}{:};
ndim         = str2double(ndim);
nno_por_elem = str2double(nno_por_elem);

switch ElemType
case 'Hexahedra'
   if nno_por_elem == 20
      idx_LaG = [1 9 2 10 3 11 4 12 13 14 15 16 5 17 6 18 7 19 8 20];
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem);
   end
case 'Tetrahedra'
   if nno_por_elem == 10
      idx_LaG = [1 5 2 6 3 7 8 9 10 4];
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem);
   end
otherwise
   fclose(fid);
   error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem);
end

fgetl(fid);  % Coordinates

xnod = zeros(10000, ndim+1);
i = 0;
while true
   linea = fgetl(fid); %  965  10  0  0
   if strcmp(linea,'end coordinates')
      break
   end;
   i = i+1;
   xnod(i,:) = sscanf(linea, '%g', [1 ndim+1]);
end
xnod = xnod(1:i, 2:ndim+1);

fgetl(fid);  % espacio en blanco
fgetl(fid);  % Elements

LaG = zeros(10000, nno_por_elem+1);
i = 0;
while true
   linea = fgetl(fid); %  965  10  0  0
   if strcmp(linea, 'end elements')
      break
   end;
   i = i+1;
   LaG(i,:) = sscanf(linea, '%d', [1 nno_por_elem+1]);
end
fclose(fid);

LaG = LaG(1:i,2:nno_por_elem+1);

LaG = LaG(:,idx_LaG);

%% bye bye!
return;
%%
