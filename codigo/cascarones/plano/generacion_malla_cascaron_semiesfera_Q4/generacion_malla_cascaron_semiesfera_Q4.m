%% Generación de la malla para elementos rectangulares de 4 nodos
clear, clc, close all

%% Se especifica el numero de elementos de la malla
nef_theta = 20; 
nef_phi   = 40;
  
%% Se definen algunas constantes
X = 1; Y = 2; Z = 3;
R = 10;

%% Eje de las theta y de las phi
theta = linspace(18, 90, nef_theta+1);
phi   = linspace(0, 360, nef_phi+1);

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
nno = (nef_theta+1)*(nef_phi+1);
nod = reshape(1:nno, nef_phi+1, nef_theta+1);

%% Se funden en una sola phi=0 y phi=360
nodos_a_remover = nod(end,:);
nod(end,:) = nod(1,:);

%% Se define la matriz xnod
xnod = [xx(:) yy(:) zz(:)];

%% Se define la matriz LaG
k = 0;
nef = nef_phi*nef_theta;
LaG = zeros(nef, 4);
for i = 1:nef_phi
   for j = 1:nef_theta
      k = k+1;
      LaG(k,:) = [    ...
         nod(i,  j)   ...  % nodo 1
         nod(i+1,j)   ...  % nodo 2
         nod(i+1,j+1) ...  % nodo 3
         nod(i,  j+1) ];   % nodo 4
   end
end

old_nno = nno;
xnod(nodos_a_remover,:) = [];
nno = size(xnod,1);
LaG = changem(LaG,1:nno,setdiff(1:old_nno,nodos_a_remover));

%% Se grafica de nuevo la malla de EFs
figure
hold on
for i = 1:nef
   plot3(xnod(LaG(i,[1:4 1]),X), xnod(LaG(i,[1:4 1]),Y), xnod(LaG(i,[1:4 1]),Z),'*-');
end
      
text(xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));

%% Se graba el resultado a disco
filename = 'semiesfera_Q4.mat';
save(filename, 'xnod', 'nno', 'LaG', 'nef');
fprintf('Resultados grabados en el archivo %s\n', filename)
