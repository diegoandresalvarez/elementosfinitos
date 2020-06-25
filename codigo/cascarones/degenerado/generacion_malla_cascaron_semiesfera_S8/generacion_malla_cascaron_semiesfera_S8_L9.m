%% Generación de la malla para elementos 
% rectangulares de 4 nodos
% serendipitos de 8 nodos  
% lagrangianos de 9 nodos

clear, clc, close all

%% Se definen algunas constantes
X = 1; Y = 2; Z = 3;
R = 10;

%% Se especifica el numero de elementos de la malla
nef_theta = 10; 
nef_phi   = 20;

%nef_theta = 3; 
%nef_phi   = 5;

%% Se define el tipo de elemento
elemento         = 'S8';  % S8 L9
elementos_planos = false;  % true false

%% Se define el nombre del archivo y se verifica si elemento está OK
if strcmp(elemento,'S8')
   if elementos_planos   
      filename = 'semiesfera_plano_S8';
   else
      filename = 'semiesfera_curvo_S8';
   end
elseif strcmp(elemento,'L9')
   if elementos_planos   
      filename = 'semiesfera_plano_L9';
   else
      filename = 'semiesfera_curvo_L9';
   end
else 
   error('Tipo de elemento debe ser S8 o L9')
end
   
%% Eje de las theta y de las phi
theta = linspace(18, 90, 2*nef_theta+1);
phi   = linspace(0, 360, 2*nef_phi+1);

%% Se crea la malla
[thth,phph] = meshgrid(theta,phi);
xx = R*sind(thth).*cosd(phph);
yy = R*sind(thth).*sind(phph);
zz = R*cosd(thth);

%% Se dibuja la malla
figure
surf(xx,yy,zz);
axis equal

%% Se numeran los nodos
nno = (2*nef_theta+1)*(2*nef_phi+1);
nod = reshape(1:nno, 2*nef_phi+1, 2*nef_theta+1);

%% Se funden en una sola phi=0 y phi=360
nodos_a_remover = nod(end,:);
nod(end,:) = nod(1,:);

%% Se define la matriz xnod
xnod = [xx(:) yy(:) zz(:)];

%% Se define la matriz LaG
k = 0;
nef = nef_phi*nef_theta;
LaG = zeros(nef, 9);
for i = 1:nef_phi
   for j = 1:nef_theta
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

old_nno = nno;
xnod(nodos_a_remover,:) = [];
nno = size(xnod,1);
LaG = changem(LaG,1:nno,setdiff(1:old_nno,nodos_a_remover));

%% Se corrigen las coordenadas xnod en el caso de cascarones formados por 
%% placas planas
if elementos_planos
   for k = 1:nef
      xnod(LaG(k,2),:) = (xnod(LaG(k,1),:) + xnod(LaG(k,3),:))/2;
      xnod(LaG(k,4),:) = (xnod(LaG(k,3),:) + xnod(LaG(k,5),:))/2;
      xnod(LaG(k,6),:) = (xnod(LaG(k,5),:) + xnod(LaG(k,7),:))/2;
      xnod(LaG(k,8),:) = (xnod(LaG(k,7),:) + xnod(LaG(k,1),:))/2;  
      xnod(LaG(k,9),:) = (xnod(LaG(k,1),:) + xnod(LaG(k,5),:))/2;       
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
text(xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));
view(3)

%% Se graba el resultado a disco
save(filename, 'xnod', 'nno', 'LaG', 'nef');

fprintf('Resultados grabados en el archivo %s.mat \n', filename)