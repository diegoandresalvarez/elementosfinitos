% Programa para generar la malla de elementos finitos de una losa de RM con 
% EF QL9

X = 1;
Y = 2;

ancho = 2;
altura = 4;
delta = 0.5;

x = 0:delta:ancho;
y = 0:delta:altura;

[xx, yy] = meshgrid(x,y);

xnod = [xx(:), yy(:)];
nno = size(xnod, 1);

figure
hold on;
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));

LaG = zeros(2*4, 9);
for i = 1:4
    LaG(i,:) = [1 10 19 20 21 12 3 2 11] + 2*(i-1);
end

for i = 5:8
    LaG(i,:) = LaG(i-4,:) + 18;
end

nef = size(LaG,1);

%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(LaG(e,:), X));
   cgy(e) = mean(xnod(LaG(e,:), Y));
   h = text(cgx(e)+0.03, cgy(e)+0.03, num2str(e)); set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
title('Malla de elementos finitos');

save malla_losa_MEN xnod LaG
disp('xnod y LaG guardados en el archivo malla_losa_MEN.mat')