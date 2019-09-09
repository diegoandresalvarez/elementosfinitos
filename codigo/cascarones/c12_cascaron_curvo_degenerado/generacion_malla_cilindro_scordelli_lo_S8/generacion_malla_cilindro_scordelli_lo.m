%% Generación de la malla para elementos serendipitos de 8 nodos y 
%% lagrangianos de 9 nodos

clear, clc, close all

%% Se especifica el numero de elementos de la malla
nef_x = 20; 
nef_t = 10;

%% Se define el tipo de elemento
elemento         = 'S8';   % S8 L9
elementos_planos = false;  % true false

%% Se define el nombre del archivo y se verifica si elemento está OK
if strcmp(elemento,'S8')
   if elementos_planos   
      filename = 'scordelli_lo_plano_S8';
   else
      filename = 'scordelli_lo_curvo_S8';
   end
elseif strcmp(elemento,'L9')
   if elementos_planos   
      filename = 'scordelli_lo_plano_L9';
   else
      filename = 'scordelli_lo_curvo_L9';
   end
else 
   error('Tipo de elemento debe ser S8 o L9')
end
   
%% Se definen algunas constantes
X = 1; Y = 2; Z = 3;
R = 25 - 0.25/2;

%% Eje de las x y de las t
x = linspace(  0, 50, 2*nef_x+1);
t = linspace(-40, 40, 2*nef_t+1);

%% Se crea la malla
[xx,tt] = meshgrid(x,t);   
yy = R*sind(tt);
zz = R*cosd(tt);

%% Se  dibuja la malla
figure
surf(xx,yy,zz);
axis equal

%% Se numeran los nodos
nno = (2*nef_x+1)*(2*nef_t+1);
nod = reshape(1:nno, 2*nef_t+1, 2*nef_x+1);

%% Se define la matriz xnod
xnod = [xx(:) yy(:) zz(:)];

%% Se define la matriz LaG
k = 0;
nef = nef_t*nef_x;
LaG = zeros(nef, 9);
for i = 1:nef_t
   for j = 1:nef_x
      k = k+1;
      LaG(k,:) = [...
         nod(2*i-1,2*j-1) ...  % nodo 1
         nod(2*i  ,2*j-1) ...  % nodo 2
         nod(2*i+1,2*j-1) ...  % nodo 3
         nod(2*i+1,2*j  ) ...  % nodo 4
         nod(2*i+1,2*j+1) ...  % nodo 5
         nod(2*i  ,2*j+1) ...  % nodo 6
         nod(2*i-1,2*j+1) ...  % nodo 7
         nod(2*i-1,2*j)   ...  % nodo 8
         nod(2*i  ,2*j) ];     % nodo 9
   end
end

%% Se corrigen las coordenadas xnod en el caso de cascarones formados por 
%% placas planas
if elementos_planos
   for k = 1:nef
      xnod(LaG(k,2),:) = (xnod(LaG(k,1),:) + xnod(LaG(k,3),:))/2;
      xnod(LaG(k,4),:) = (xnod(LaG(k,3),:) + xnod(LaG(k,5),:))/2;
      xnod(LaG(k,6),:) = (xnod(LaG(k,5),:) + xnod(LaG(k,7),:))/2;
      xnod(LaG(k,8),:) = (xnod(LaG(k,7),:) + xnod(LaG(k,1),:))/2;  
   end
end

%% Si se tiene un elemento serendipito de 8 nodos
if strcmp(elemento, 'S8')
   % Se identifican los nodos a eliminar y se identifican sus indices
   nodos_elim = unique(LaG(:,9))';
   old_idx = setdiff(1:nno,nodos_elim);
   
   % Se elimina la columna 9 de LaG
   LaG(:,9) = [];
   
   % Se eliminan los nodos redundantes de la matriz LaG
   xnod(nodos_elim, :) = [];
   nno = size(xnod,1);
   
   % Se renumera la matriz LaG
   LaG = changem(LaG, 1:nno, old_idx);
end

%% Se grafica de nuevo la malla de EFs
figure
hold on
switch elemento
   case 'S8'
      for i = 1:nef
         plot3(xnod(LaG(i,[1:8 1]),X), xnod(LaG(i,[1:8 1]),Y), xnod(LaG(i,[1:8 1]),Z),'*-');
      end
   case 'L9'
      for i = 1:nef
         plot3(xnod(LaG(i,[1:8 1]),X), xnod(LaG(i,[1:8 1]),Y), xnod(LaG(i,[1:8 1]),Z),'*-');
         plot3(xnod(LaG(i,9),X), xnod(LaG(i,9),Y), xnod(LaG(i,9),Z),'*');         
      end
end   
view(3)

%% Se graba el resultado a disco
save(filename, 'xnod', 'nno', 'LaG', 'nef');

fprintf('Resultados grabados en el archivo %s.mat \n', filename)