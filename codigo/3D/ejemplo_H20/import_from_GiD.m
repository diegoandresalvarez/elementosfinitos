function [xnod,LaG,material] = import_from_GiD(filename)
% Esta rutina importa los archivos de malla generados por GiD
%
% Uso:
% [xnod,LaG]          = import_from_GiD(filename)
% [xnod,LaG,material] = import_from_GiD(filename)
%
% VER EN LA AYUDA DE GID:
% Postprocess mesh format: ProjectName.post.msh
% Postprocess results format: ProjectName.post.res
% http://www-opale.inrialpes.fr/Aerochina/info/en/html-version/gid_17.html

%% SE ABRE EL ARCHIVO
fid = fopen([filename '.msh'], 'r');

%% SE LEEN LOS COMENTARIOS Y LINEAS BLANCAS INICIALES
linea = leer_linea(fid);

%% LINEA MESH
%       mesh dimension 3 elemtype hexahedra nnode 20
expr = 'MESH (\w*)DIMENSION (\d+) ELEMTYPE (\w+) NNODE (\d+)';

tmp = regexp(linea, expr, 'tokens');
[~, ndim, ElemType, nno_por_elem] = tmp{1}{:};
ndim = str2double(ndim);
if ndim ~= 2 && ndim ~= 3
    error('DIMENSION debe ser 2 o 3');
end

nno_por_elem = str2double(nno_por_elem);

switch ElemType
case 'POINT'
   if nno_por_elem == 1
      idx_LaG = 1;
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end
case 'LINEAR'
   if nno_por_elem == 2 
      idx_LaG = [ 1 2 ];
   elseif nno_por_elem == 3
      idx_LaG = [ 1 3 2 ];       
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end   
case 'TRIANGLE'
   if nno_por_elem == 3
      idx_LaG = [ 1 2 3 ];             
   elseif  nno_por_elem == 6
      idx_LaG = [ 1 4 2 5 3 6 ];                   
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end      
case 'QUADRILATERAL'
   if nno_por_elem == 4
      idx_LaG = [ 1 2 3 4 ];             
   elseif  nno_por_elem == 8
      idx_LaG = [ 1 5 2 6 3 7 4 8 ];
   elseif  nno_por_elem == 9
      idx_LaG = [ 1 5 2 6 3 7 4 8 9 ];
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end
case 'TETRAHEDRA'
   if nno_por_elem == 4
      idx_LaG = [ 1 2 3 4 ];
   elseif nno_por_elem == 10
      idx_LaG = [ 1 5 2 6 3 7 8 9 10 4 ];
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end   
case 'HEXAHEDRA'
   if nno_por_elem == 8
      idx_LaG = [ 1 2 3 4 5 6 7 8 ];
   elseif nno_por_elem == 20
      idx_LaG = [ 1 9 2 10 3 11 4 12 13 14 15 16 5 17 6 18 7 19 8 20 ];
   elseif nno_por_elem == 27
      idx_LaG = [ 1 9 2 10 3 11 4 12 13 22 14 23 15 24 16 25 5 17 6 18 7 19 8 20 ];
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end
case 'PRISM'
   if nno_por_elem == 4
      idx_LaG = NaN;
   elseif nno_por_elem == 15
      idx_LaG = NaN;
   else
      fclose(fid);
      error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem)
   end     
%case {'pyramid', 'sphere', 'circle'}
otherwise
   fclose(fid);
   error('Elemento %s de %d nodos no soportado', ElemType, nno_por_elem);
end

linea = leer_linea(fid);

%% COORDINATES .. END COORDINATES
if not(strcmpi(linea, 'COORDINATES'))
    error('Se espera el bloque COORDINATES .. END COORDINATES')
end

xnod = zeros(10000, ndim+1);
i = 0;
while true
   linea = leer_linea(fid); %  965  10  0  0
   if strcmpi(linea,'END COORDINATES')
      break
   end;
   i = i+1;
   xnod(i,:) = sscanf(linea, '%g', [1 ndim+1]);
end

% se recorta la matriz xnod
xnod = xnod(1:i,:);

% se escriben las filas de acuerdo con los subindices
idx = xnod(:,1);
xnod = xnod(idx, 2:ndim+1);

linea = leer_linea(fid);

%% ELEMENTS .. END ELEMENTS
if not(strcmpi(linea, 'ELEMENTS'))
    error('Se espera el bloque ELEMENTS .. END ELEMENTS')
end

% se lee la primera linea
linea = leer_linea(fid); 
% #el_num node_1 node_2 node_3 material (optional)
linea1 = sscanf(linea, '%d', [1 nno_por_elem+2]);
if nno_por_elem ~= length(linea1)-1 && nno_por_elem ~= length(linea1)-2
    error('El numero de columnas de ELEMENTS no coincide con el especificado')
end

% ya se sabe el numero de columnas con el cual crear LaG
LaG = zeros(10000, length(linea1));
LaG(1,:) = linea1;
i = 1;
while true
   % #el_num node_1 node_2 node_3 material (optional)
   %  965  10  0  0 1
   linea = leer_linea(fid); 
   if strcmpi(linea, 'END ELEMENTS')
      break
   end;
   i = i+1;
   LaG(i,:) = sscanf(linea, '%d', [1 nno_por_elem+2]);
end
fclose(fid);

% se recorta la matriz LaG
LaG = LaG(1:i,:);

% se escriben las filas de acuerdo con los subindices
idx = LaG(:,1);
if size(LaG,2) == nno_por_elem+2
    material(idx,:) = LaG(:,nno_por_elem+2);
elseif size(LaG,2) == nno_por_elem+1
    material = NaN;
else
    error('El numero de columnas de ELEMENTS no coincide con el especificado')
end
LaG = LaG(idx,2:nno_por_elem+1);

% se reorganizan las columnas con mi numeracion local de los nodos
LaG = LaG(:,idx_LaG);

%% bye bye!
end

%%
function linea = leer_linea(fid)
% Lea una línea ignorando las líneas en blanco iniciales y los comentarios
    while true
        linea = strtrim(upper(fgetl(fid)));
        if not(isempty(linea)) && linea(1) ~= '#'
            break
        end
    end
end
