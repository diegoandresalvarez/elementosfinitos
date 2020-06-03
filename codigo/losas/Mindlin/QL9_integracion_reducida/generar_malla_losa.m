%% Creacion de la malla de EFs de una losa de Mindlin con el EF QL9

%% Seleccione la malla a generar
% malla con delta = 0.1 para el ejemplo general
%malla = 'malla_losa';     

% malla con delta = 0.5 para graficar los modos de energia nula
malla = 'malla_losa_MEN';  

%% Se definen algunas constantes para facilitar la lectura del codigo
X = 1; Y = 2;

%% Dimensiones de la losa y tamanio de cada EF
ancho  = 2;
altura = 4;
switch malla
    case 'malla_losa_MEN'
        delta  = 0.5;
    case 'malla_losa'
        delta  = 0.1;
end

%% Matriz de coordenadas de los nodos xnod
x = 0:delta:ancho;
y = 0:delta:altura;

[xx, yy] = meshgrid(x,y);

xnod = [xx(:), yy(:)];
nno = size(xnod, 1);   % numero de nodos

%% Matriz de conectividades LaG
switch malla
    case 'malla_losa_MEN'
        LaG = zeros(2*4, 9);
        for i = 1:4
            LaG(i,:) = [1 10 19 20 21 12 3 2 11] + 2*(i-1);
        end

        for i = 5:8
            LaG(i,:) = LaG(i-4,:) + 18;
        end
    case 'malla_losa'
        LaG = zeros(10*20, 9);
        for i = 1:20
            LaG(i,:) = [1 42 83 84 85 44 3 2 43] + 2*(i-1);
        end

        k = 21;
        for j = 2:10
            for i = 1:20
                LaG(k,:) = LaG(i,:) + 82*(j-1);
                k = k+1;
            end
        end        
end

nef = size(LaG,1);   % numero de EFs

%% Se calcula el centro de gravedad de cada EF
cgx = zeros(1,nef); 
cgy = zeros(1,nef);
for e = 1:nef
   cgx(e) = mean(xnod(LaG(e,:), X));
   cgy(e) = mean(xnod(LaG(e,:), Y));
end

%% Se dibuja la malla de elementos finitos
figure
hold on
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y)); 
   text(cgx(e)+0.03, cgy(e)+0.03, num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'rx');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
axis([-delta, ancho+delta, -delta, altura+delta])
title('Malla de una losa con EFs QL9');

%% Se guarda la malla en disco
disp(['xnod y LaG guardados en el archivo "' malla '.mat".'])
% save(malla, 'xnod', 'LaG');