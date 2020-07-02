%% TECHO DE SCORDELIS-LO
%% Generacion de la malla para elementos rectangulares de 4 nodos


clear, clc, close all

%% Se especifica el numero de elementos de la malla
nef_x = 30; 
nef_t = 10;

%% Se definen algunas constantes
X = 1; Y = 2; Z = 3;
R = 25;

%% Eje de las x y de las t
x = linspace(  0, 50, nef_x+1);
t = linspace(-40, 40, nef_t+1);

%% Se crea la malla
[xx,tt] = meshgrid(x,t);   
yy = R*sind(tt);
zz = R*cosd(tt);

%% Se  dibuja la malla
figure
surf(xx,yy,zz);
axis equal

%% Se numeran los nodos
nno = (nef_x+1)*(nef_t+1);
nod = reshape(1:nno, nef_t+1, nef_x+1);

%% Se define la matriz xnod
xnod = [xx(:) yy(:) zz(:)];

%% Se define la matriz LaG
k = 0;
nef = nef_t*nef_x;
LaG = zeros(nef, 4);
for i = 1:nef_t
   for j = 1:nef_x
      k = k+1;
      LaG(k,:) = [    ...
         nod(i,  j)   ...  % nodo 1
         nod(i+1,j)   ...  % nodo 2
         nod(i+1,j+1) ...  % nodo 3
         nod(i,  j+1) ];   % nodo 4
   end
end

%% Se grafica de nuevo la malla de EFs
figure
hold on
for i = 1:nef
   plot3(xnod(LaG(i,[1:4 1]),X), xnod(LaG(i,[1:4 1]),Y), xnod(LaG(i,[1:4 1]),Z),'*-');
end
      
text(xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));
view(3)

%% Se graba el resultado a disco
filename = sprintf('scordelli_lo_Q4_malla_%d_%d', nef_x, nef_t);
save(filename, 'xnod', 'nno', 'LaG', 'nef');
fprintf('Resultados grabados en el archivo %s.mat \n', filename);
