clear, clc, close all % borro la memoria, la pantalla y las figuras

%% ------------------------------------------------------------------------
%% NOTA: este codigo SOLO es apropiado para TENSION PLANA
%% ------------------------------------------------------------------------

%% definicion del problema
% Calcule los desplazamientos y las reacciones en los empotramiento, las
% deformaciones y los esfuerzos de la estructura en TENSION PLANA mostrada 
% en la figura adjunta

%% defino las variables/constantes
X    = 1;           % un par de constantes que ayudaran en la 
Y    = 2;           % lectura del codigo
Ee   = 200e9;       % modulo de elasticidad del solido (Pa) = 200GPa
nue  = 0.30;        % coeficiente de Poisson
te   = 0.10;        % espesor del solido (m)
rhoe = 7850;        % densidad (kg/m^3)
g    = 9.81;        % aceleracion de la gravedad (m/s^2)

% Malla ejemplo
%malla_ejemplo

% Malla_refinada (malla elaborada por David Felipe Cano Perdomo)
malla_refinada


%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1 2 3 1]),X), xnod(LaG(e,[1 2 3 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X) + xnod(LaG(e,3),X))/3;
   cgy(e) = (xnod(LaG(e,1),Y) + xnod(LaG(e,2),Y) + xnod(LaG(e,3),Y))/3;   
   h = text(cgx(e), cgy(e), num2str(e)); set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
B = cell(nef,1);       % contenedor para las matrices de deformacion

% matriz constitutiva del elemento para TENSION PLANA
De = [ Ee/(1-nue^2)     Ee*nue/(1-nue^2)  0
       Ee*nue/(1-nue^2) Ee/(1-nue^2)      0
       0                0                 Ee/(2*(1+nue)) ];

for e = 1:nef      % ciclo sobre todos los elementos finitos
   % Calculo de la matriz de rigidez del elemento e
   x1 = xnod(LaG(e,1),X);              y1 = xnod(LaG(e,1),Y);
   x2 = xnod(LaG(e,2),X);              y2 = xnod(LaG(e,2),Y);
   x3 = xnod(LaG(e,3),X);              y3 = xnod(LaG(e,3),Y);
   
   Ae = 0.5*det([ 1 x1 y1      %Area del EF e
                  1 x2 y2
                  1 x3 y3]);               
   if Ae <= 0
      error('Revise las coordenadas locales del EF %d.\n', e);
   end
   
   % Calculo de la matriz de deformaciones B.
   a1 = x2*y3 - x3*y2;        b1 = y2-y3;        c1 = x3-x2;
   a2 = x3*y1 - x1*y3;        b2 = y3-y1;        c2 = x1-x3;
   a3 = x1*y2 - x2*y1;        b3 = y1-y2;        c3 = x2-x1;
   
   B{e} = (1/(2*Ae))*[ b1    0      b2    0      b3    0 
                        0   c1       0   c2       0   c3
                       c1   b1      c2   b2      c3   b3 ];
   
   Ke = B{e}'*De*B{e}*te*Ae;
   
   % Calculo del vector de fuerzas nodales equivalentes del elemento e
   % Fuerzas masicas (peso propio)
   fbe = [0; -rhoe*g; 0; -rhoe*g; 0; -rhoe*g] * Ae*te/3;
        
   fe = fbe; % vector de fuerzas nodales equivalentes
   
   % Ensamblo las contribuciones a las matrices globales
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) ];
   K(idx,idx) = K(idx,idx) + Ke;
   f(idx,:)   = f(idx,:)   + fe;
end;

%% Relacion de las cargas superficiales (vector ft)
ft = sparse(ngdl,1); % fuerzas nodales equivalentes de cargas superficiales
for i = 1:nlcd
   e     = carga_distr(i,1);
   lado  = carga_distr(i,2);
   carga = carga_distr(i,3:6);
   fte = t2ft_T3(xnod(LaG(e,[1 2 3]),[X Y]), lado, carga, te);
   
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) ];
   ft(idx,:) = ft(idx,:) + fte;
end

% Agrego al vector de fuerzas nodales equivalentes las fuerzas
% superficiales calculadas
f = f + ft;

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos
c = restric(:,1);   d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |  % recuerde que qc=0 (siempre)
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = restric(:,2);   % desplaz. conocidos (para gdl(1,X) gdl(1,Y) gdl(5,Y))

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% imprimo los resultados
format short g
disp('Nodo   Despl_x (m)   Despl_y (m) = ');     [1:nno; reshape(a,2,nno)]'
disp('Nodo Fuerzas nodales equiv. X, Y (N) = '); [1:nno; reshape(f,2,nno)]'
disp('Nodo Fuerzas nodales equil. X, Y (N) = '); [1:nno; reshape(q,2,nno)]'

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,2,nno)';
escala = 20000; % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   h1 = line(xnod(LaG(e,[1 2 3 1]),X), xnod(LaG(e,[1 2 3 1]),Y)); %original
   set(h1, 'Color', [0 0 1]); % color expresado en notacion RBG entre 0 y 1
   h2 = line(xdef(LaG(e,[1 2 3 1]),X), xdef(LaG(e,[1 2 3 1]),Y)); %deformada
   set(h2, 'Color', [1 0 0]);
end
axis equal tight
legend('Posicion original','Posicion deformada','Location', 'SouthOutside')
title(sprintf('Deformada escalada %d veces',escala));

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = zeros(3,nef);
esf = zeros(3,nef);
for e = 1:nef
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) ];     
   ae = a(idx);            % desplazamientos de los gdl del elemento e
   def(:,e) = B{e}*ae;     % Calculo las deformaciones
   esf(:,e) = De*def(:,e); % Calculo los esfuerzos
end;
sx = esf(1,:);  sy = esf(2,:);  txy = esf(3,:);
ex = def(1,:);  ey = def(2,:);  gxy = def(3,:); 
ez  = -(nue/Ee)*(sx+sy); % Se calculan las deformacion ez en tension plana

%% imprimo y grafico las deformaciones
disp('Deformaciones: (EF,ex,ey,ez,gxy) = '); [1:nef; ex; ey; ez; gxy]'
figure
subplot(4,1,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ex(e))
end;
ylabel('\epsilon_x','FontSize',26); axis equal tight; colorbar; 

subplot(4,1,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ey(e))
end;
ylabel('\epsilon_y','FontSize',26); axis equal tight; colorbar;

subplot(4,1,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ez(e))
end;
ylabel('\epsilon_z','FontSize',26); axis equal tight; colorbar;

subplot(4,1,4); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),gxy(e))
end;
ylabel('\gamma_{xy}','FontSize',26); axis equal tight; colorbar;

%% imprimo y grafico los esfuerzos
disp('Esfuerzos (Pa):  (EF,sx,sy,txy) = '); [1:nef; sx; sy; txy]'
figure
subplot(3,1,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sx(e))
end;
ylabel('\sigma_x (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(3,1,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sy(e))
end;
ylabel('\sigma_y (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(3,1,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),txy(e))
end;
ylabel('\tau_{xy} (Pa)','FontSize',26); axis equal tight; colorbar;


%% Se calculan y grafican para cada elemento los esfuerzos principales y
%% sus direcciones
% NOTA: esto solo es valido para el caso de TENSION PLANA). 
% En caso de DEFORMACION PLANA se deben calcular los valores y vectores 
% propios de la matriz de tensiones de Cauchy
%   [dirppales{e}, esfppales{e}] = eig([sx  txy 0    % matriz de esfuerzos
%                                       txy sy  0    % de Cauchy
%                                       0   0   0]);

s1   = (sx+sy)/2 + sqrt(((sx-sy)/2).^2+txy.^2); % esfuerzo normal maximo
s2   = (sx+sy)/2 - sqrt(((sx-sy)/2).^2+txy.^2); % esfuerzo normal minimo
tmax = (s1-s2)/2;                               % esfuerzo cortante maximo
ang  = 0.5*atan2(2*txy, sx-sy); % angulo de inclinacion de s1

%% Calculo de los esfuerzos de von Mises
s3 = zeros(size(s1));
sv = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2);

%% imprimo los resultados
disp('Elemento,s1(Pa),s2(Pa),tmax(Pa),angulo(rad) = '); [1:nef; s1; s2; tmax; ang]'
disp('Elemento,Esfuerzos de von Mises (Pa) = '); [1:nef; sv]'

esc = 2; % escala para graficar las flechas
figure
subplot(3,1,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s1(e))
end;
% Grafique lineas que indiquen direcciones principales de sigma_1
quiver(cgx,cgy,...     %  En el punto (cgx,cgy) grafique una flecha (linea)
   s1.*cos(ang),s1.*sin(ang),... % indicando la direccion principal de sigma_1
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(cgx,cgy,...               % La misma flecha ahora en la otra direccion,
   s1.*cos(ang+pi),s1.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\sigma_1 (Pa)','FontSize',26); colorbar

subplot(3,1,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s2(e))
end;
% Grafique lineas que indiquen direcciones principales de sigma_2
quiver(cgx,cgy,...                         % flecha indicando la direccion 
   s2.*cos(ang+pi/2),s2.*sin(ang+pi/2),... % principal de sigma_2
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, s2.*cos(ang-pi/2),s2.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\sigma_2 (Pa)','FontSize',26); colorbar

subplot(3,1,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),tmax(e))
end;
% Grafique lineas que indiquen direcciones principales de tau_max,
quiver(cgx,cgy, tmax.*cos(ang+pi/4),tmax.*sin(ang+pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang-pi/4),tmax.*sin(ang-pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang+3*pi/4),tmax.*sin(ang+3*pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(cgx,cgy, tmax.*cos(ang-3*pi/4),tmax.*sin(ang-3*pi/4),'k',...
         'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal; axis([-0.1 0.9 -0.1 0.4]);
ylabel('\tau_{max} (Pa)','FontSize',26); colorbar

figure; hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sv(e))
end;
ylabel('\sigma_v (Pa)','FontSize',26); axis equal tight; colorbar;
title('Esfuerzos de von Mises (Pa)')

%%
return; % bye, bye!
