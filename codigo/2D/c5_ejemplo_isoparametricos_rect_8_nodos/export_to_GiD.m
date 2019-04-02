function export_to_GiD(filename,xnod,LaG,a,q,stress)
% Esta rutina genera los archivos requeridos para exportar los resultados
% de MATLAB a GiD.
%
% Uso:
% export_to_GiD(filename,xnod,LaG,a,q,stress)
%
% donde:
%
%
%
%

% VER EN LA AYUDA DE GID:
% Postprocess mesh format: ProjectName.post.msh
% Postprocess results format: ProjectName.post.res
% http://www-opale.inrialpes.fr/Aerochina/info/en/html-version/gid_17.html

%% Inicializar variables
nno          = size(xnod,1); % numero de nodos (numero de filas de xnod)
ndim         = size(xnod,2); % numero de dimensiones del problema
nef          = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nno_por_elem = size(LaG,2);  % numero de nodos por elemento

if ndim == 2
   switch nno_por_elem
      case {3,6}
         eletyp = 'Triangle';
      case 4
         eletyp = 'Quadrilateral';
      case 8
         eletyp = 'Quadrilateral';
         LaG = LaG(:, [1 3 5 7 2 4 6 8]);
      case 9
         eletyp = 'Quadrilateral';
         LaG = LaG(:, [1 3 5 7 2 4 6 8 9]);
      otherwise
         error('Tipo de elemento 2D no soportado');
   end

   if isa(stress,'cell')
      npg = size(esf,2);  % numero de puntos de Gauss
      switch nno_por_elem
         case {3,6}   % 'Triangle';
            if ~ismember(npg,[1 3 6])
               error('Numero de puntos de Gauss invalido');
            end
         case {4,8,9} % 'Quadrilateral';
            if ~ismember(npg,[1 4 9])
               error('Numero de puntos de Gauss invalido');
            end
      end
   end
elseif ndim == 3
   switch nno_por_elem
      case 4
         eletyp = 'Tetrahedron'; % ?? Tetraedra?
      case 10
         eletyp = 'Tetrahedron';
         %LaG = LaG(:, [FALTA]);
      case 8
         eletyp = 'Hexahedron';  % ?? Hexahedra?
      case 20
         eletyp = 'Hexahedra';
         LaG(:,[1 9 2 10 3 11 4 12 13 14 15 16 5 17 6 18 7 19 8 20]) = LaG;
      case 27
         eletyp = 'Hexahedron';  % ?? Hexahedra?
         %LaG = LaG(:, [FALTA]);
      otherwise
         error('Tipo de elemento 3D no soportado');
   end
   
   if isa(stress,'cell')
      npg = size(stress,2);  % numero de puntos de Gauss

      switch nno_por_elem
         case {4,10}    % 'Tetrahedron'
            if ~ismember(npg,[1 4 10])
               error('Numero de puntos de Gauss invalido');
            end
         case {8,20,27} % 'Hexahedron';  % ?? Hexahedra?
            if ~ismember(npg,[1 8 27])
               error('Numero de puntos de Gauss invalido');
            end
      end
   end
else
   error('Numero de dimensiones invalida');
end

%% ***** ***** ***** Escribiendo la malla ***** ***** *****
fid = fopen([filename '.msh'],'w');

%% Escribiendo encabezado del archivo:
fprintf(fid,'### \n');
fprintf(fid,'# MAT_FEM  V.1.0 \n');
fprintf(fid,'# \n');
fprintf(fid,'MESH dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n', ndim, eletyp, nno_por_elem);

%% Escribiendo coordenadas de los nodos:
fprintf(fid,'coordinates \n');
for i = 1:nno
   fprintf(fid, ['%5d ' repmat('%12.5f ',1,ndim) '\n'], i, xnod(i,:));
end
fprintf(fid,'end coordinates \n \n');

%% Escribiendo elementos:
fprintf(fid,'elements \n');
for i = 1:nef
   % El uno final es "material number"
   fprintf(fid,['%5d ' repmat('%5d ',1,nno_por_elem) '     1 \n'],i,LaG(i,:));
end
fprintf(fid,'end elements \n \n');
fclose(fid);

%% ***** ***** ***** Escribiendo los resultados ***** ***** *****
fid = fopen([filename '.res'], 'w');

%% Escribiendo encabezado de la malla:
fprintf(fid,'Gid Post Results File 1.0 \n');
fprintf(fid,'### \n');
fprintf(fid,'# MAT_FEM  V.1.0 \n');
fprintf(fid,'# \n\n');

%% Escribiendo los desplazamientos:
% Result "result name" "analysis name" step_value my_result_type my_location "location name"
%  Result:          a keyword that should be written as it is, case doesn't matter.
%  "result name":   a name for the Result, which will be used for menus.
%  "analysis name": the name of the analysis of this Result, which will be used for menus.
%  step_value:      the value of the step inside the analysis "analysis name".
%  my_type:         type of the Result, should be one of Scalar, Vector,
%                   Matrix, PlainDeformationMatrix, MainMatrix, LocalAxes.
%  my_location:     where is the Result located, should be one of OnNodes,
%                   OnGaussPoints. If the Result is OnGaussPoints a
%                   "location name" should be entered.
%  "location name": name of the Gauss Points on which the Result is
%                   defined
fprintf(fid,'Result "Displacements" "Load Analysis"  1  Vector OnNodes \n');
fprintf(fid,'ComponentNames "X-Displ", "Y-Displ", "Z-Displ" \n');
fprintf(fid,'Values \n');
switch ndim
   case 2
      for i = 1:nno
         fprintf(fid,'%5d %13.5g %13.5g %13.5g \n', i,a(i*2-1),a(i*2), 0);
      end
   case 3
      for i = 1:nno
         fprintf(fid,'%5d %13.5g %13.5g %13.5g \n', i,a(i*3-2),a(i*3-1),a(i*3));
      end
end
fprintf(fid,'End Values \n');
fprintf(fid,'# \n\n');

%% Escribiendo las fuerzas nodales de equilibrio:
fprintf(fid,'Result "Reaction" "Load Analysis"  1  Vector OnNodes \n');
fprintf(fid,'ComponentNames "Rx", "Ry", "Rz" \n');
fprintf(fid,'Values \n');
switch ndim
   case 2
      for i = 1:nno
         fprintf(fid,'%5d %12.5g %12.5g %12.5g \n', i, q(i*2-1), q(i*2), 0);
      end
   case 3
      for i = 1:nno
         fprintf(fid,'%5d %12.5g %12.5g %12.5g \n', i,q(i*3-2),q(i*3-1),q(i*3));
      end
end
fprintf(fid,'End Values \n');
fprintf(fid,'# \n\n');

if isa(stress,'double')
   %% Escribiendo los esfuerzos:
   fprintf(fid,'Result "Stresses" "Load Analysis"  1  Matrix OnNodes \n');
   fprintf(fid,'ComponentNames "sx", "sy", "sz", "txy", "txz", "tyz" \n');
   fprintf(fid,'Values \n');
   for i = 1:nno
      %            i   sx     sy     sz     txy    txz    tyz
      fprintf(fid,'%5d %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g \n',i,stress(i,:));
   end
   fprintf(fid,'End Values \n');
elseif isa(stress,'cell')
   %% Escribiendo los puntos de Gauss:
   fprintf(fid,'GaussPoints "gauss_points"    Elemtype %s \n', eletyp);
   fprintf(fid,'Number of Gauss Points: %2.0f \n',npg);
   fprintf(fid,'Natural Coordinates: Internal \n');
   fprintf(fid,'End GaussPoints \n');
   fprintf(fid,'# \n\n');
   
   %% Escribiendo los esfuerzos:
   fprintf(fid,'Result "Stresses" "Load Analysis"  1  Matrix OnGaussPoints "gauss_points" \n');
   fprintf(fid,'ComponentNames "sx", "sy", "sz", "txy", "txz", "tyz" \n');
   fprintf(fid,'Values \n');
   switch eletyp
      case 'Quadrilateral'
         if npg == 4
            pos_nodos = [ ...
               1 1               %  IV   III
               2 1               %  I    II
               2 2
               1 2 ];
         else
            error('No se ha definido la posicion relativa de los puntos de Gauss');
         end
      case 'Hexahedra'
         if npg == 8
            pos_nodos = [ ...  % GiD en nivel abajo     GiD en nivel arriba
               1               %  IV   III              VIII   VII
               5               %  I    II               V      VI
               7
               3
               2               % mi programa            mi programa
               6               % III   VII              IV     VIII
               8               % I     V                II     VI
               4 ];
         else
            error('No se ha definido la posicion relativa de los puntos de Gauss');
         end
      otherwise
         error('No se ha definido la posicion relativa de los puntos de Gauss');         
   end

   for i = 1:nef
      for j = 1:npg
         if j == 1
            %            i
            fprintf(fid,'%5d',i);
         else
            fprintf(fid,'     ');
         end;
         
         if ndim == 2
         %             sx     sy     sz     txy    txz    tyz
         fprintf(fid,' %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g \n', ...
            stress{i,pos_nodos(j,1),pos_nodos(j,2)}(1), ...
            stress{i,pos_nodos(j,1),pos_nodos(j,2)}(2), ...
            0,                                          ...
            stress{i,pos_nodos(j,1),pos_nodos(j,2)}(3), ...
            0,                                          ...
            0);
         elseif ndim == 3
         %             sx     sy     sz     txy    txz    tyz            
         fprintf(fid,' %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g \n', ...
            stress{i,pos_nodos(j)}(1), ...
            stress{i,pos_nodos(j)}(2), ...
            stress{i,pos_nodos(j)}(3), ...
            stress{i,pos_nodos(j)}(4), ...
            stress{i,pos_nodos(j)}(5), ...
            stress{i,pos_nodos(j)}(6));
         end
      end
   end
   fprintf(fid,'End Values \n');
end
fclose(fid);

%% bye bye!
return;
%%