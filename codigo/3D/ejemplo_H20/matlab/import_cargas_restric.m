function [cargas,restric] = import_cargas_restric(filename)
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
fid = fopen([filename '.cr'], 'r');

%% SE LEEN LOS COMENTARIOS Y LINEAS BLANCAS INICIALES
linea = leer_linea(fid);

%% COORDINATES .. END COORDINATES
if not(strcmpi(linea, 'CARGAS'))
    error('Se espera el bloque CARGAS .. END CARGAS')
end

cargas = zeros(10000, 3);
i = 0;
while true
   linea = leer_linea(fid);
   if strcmpi(linea,'END CARGAS')
      break
   end
   i = i+1;
   cargas(i,:) = sscanf(linea, '%g', [1 3]);
end

% se recorta la matriz xnod
cargas = cargas(1:i,:);

linea = leer_linea(fid);

%% ELEMENTS .. END ELEMENTS
if not(strcmpi(linea, 'RESTRICCIONES'))
    error('Se espera el bloque RESTRICCIONES .. END RESTRICCIONES')
end

restric = zeros(10000, 3);
i = 0;
while true
   linea = leer_linea(fid); 
   if strcmpi(linea, 'END RESTRICCIONES')
      break
   end
   i = i+1;
   restric(i,:) = sscanf(linea, '%g', [1 3]);
end
fclose(fid);

% se recorta la matriz LaG
restric = restric(1:i,:);

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
