% Programa para generar la malla de elementos finitos

X = 1;
Y = 2;

ancho = 2;
altura = 4;
delta = 0.1;

x = 0:delta:ancho;
y = 0:delta:altura;

[xx, yy] = meshgrid(x,y);

xnod = [xx(:), yy(:)];
nno = size(xnod, 1);

% figure
% hold on;
% plot(xnod(:,X), xnod(:,Y), 'r*');
% text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));

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

save malla_losa xnod LaG
disp('xnod y LaG guardados en el archivo malla_losa.mat')