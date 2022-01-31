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
malla_refinada_v1

%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(nef,1); cgy = zeros(nef,1); % almacena el centro de gravedad
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
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
B   = cell(nef,1);       % contenedor para las matrices de deformacion
idx = cell(nef,1);       % gdl de cada EF (para el ensamblaje)

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
   
   % Calculo de la matriz de deformaciones B del EF e
   b1 = y2-y3;        c1 = x3-x2;
   b2 = y3-y1;        c2 = x1-x3;
   b3 = y1-y2;        c3 = x2-x1;
   
   B{e} = (1/(2*Ae))*[ b1    0      b2    0      b3    0 
                        0   c1       0   c2       0   c3
                       c1   b1      c2   b2      c3   b3 ];
   
   % Calculo de la matriz de rigidez del EF e
   Ke = te*B{e}'*De*B{e}*Ae;
   
   % Calculo del vector de f.n.e. de fuerzas masicas del EF e (peso propio)
   fbe = -rhoe*g*Ae*te*[0; 1; 0; 1; 0; 1]/3;

   % Ensamblo las contribuciones a las matrices globales
   idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) ];
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fbe;
end

%% Calculo del vector de f.n.e. de fuerzas superficiales del EF e
for i = 1:nlcd
   e     = carga_distr(i,1);
   lado  = carga_distr(i,2);
   carga = carga_distr(i,3:6);
   fte = t2ft_T3(xnod(LaG(e,[1 2 3]),[X Y]), lado, carga, te);

   f(idx{e},:) = f(idx{e},:) + fte;
end

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
ac = restric(:,2);   % desplazamientos conocidos

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
escala = 20000;             % factor de escalamiento de la deformada
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
   ae = a(idx{e});         % desplazamientos de los gdl del elemento e
   def(:,e) = B{e}*ae;     % Calculo las deformaciones
   esf(:,e) = De*def(:,e); % Calculo los esfuerzos
end
sx = esf(1,:)';  sy = esf(2,:)';  txy = esf(3,:)';
ex = def(1,:)';  ey = def(2,:)';  gxy = def(3,:)'; 
ez  = -(nue/Ee)*(sx+sy); % Se calculan las deformacion ez en tension plana

%% imprimo y grafico las deformaciones
disp('Deformaciones: (EF,ex,ey,ez,gxy) = '); [(1:nef)' ex ey ez gxy]
figure
subplot(2,2,1); plot_def_esf(xnod, LaG, ex,  '\epsilon_x');
subplot(2,2,2); plot_def_esf(xnod, LaG, ey,  '\epsilon_y');
subplot(2,2,3); plot_def_esf(xnod, LaG, ez,  '\epsilon_z');
subplot(2,2,4); plot_def_esf(xnod, LaG, gxy, '\gamma_{xy} [rad]');

%% imprimo y grafico los esfuerzos
disp('Esfuerzos (Pa):  (EF,sx,sy,txy) = '); [(1:nef)' sx sy txy]
figure
subplot(3,1,1); plot_def_esf(xnod, LaG, sx,  '\sigma_x [Pa]');
subplot(3,1,2); plot_def_esf(xnod, LaG, sy,  '\sigma_y [Pa]');
subplot(3,1,3); plot_def_esf(xnod, LaG, txy, '\tau_{xy} [Pa]');

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
disp('Elemento,s1(Pa),s2(Pa),tmax(Pa),angulo(rad) = '); [(1:nef)' s1 s2 tmax ang]
disp('Elemento,Esfuerzos de von Mises (Pa) = '); [(1:nef)' sv]

figure
subplot(2,2,1); plot_def_esf(xnod, LaG, s1,   '(\sigma_1)_{xy} [Pa]', cgx, cgy, { ang })
subplot(2,2,2); plot_def_esf(xnod, LaG, s2,   '(\sigma_2)_{xy} [Pa]', cgx, cgy, { ang+pi/2 })
subplot(2,2,3); plot_def_esf(xnod, LaG, tmax, '\tau_{max} [Pa]',      cgx, cgy, { ang+pi/4, ang-pi/4 })
subplot(2,2,4); plot_def_esf(xnod, LaG, sv,   'Esfuerzos de von Mises [Pa]');

%%
return; % bye, bye!

function plot_def_esf(xnod, LaG, variable, texto, cgx, cgy, angulos)
    X = 1; Y = 2;
    hold on; 
    colorbar;
    
    nef = size(LaG, 1);    
    for e = 1:nef  
       fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), variable(e));
    end
    axis equal tight
    colormap jet
    title(texto, 'FontSize',20);
   
    esc = 0.5;
    if nargin == 7
        norma = 1; % = variable % si se quiere proporcional
        for i = 1:length(angulos)
            % se indica la flecha de la direccion principal
            quiver(cgx, cgy,...             
                norma.*cos(angulos{i}), norma.*sin(angulos{i}),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off',... % una flecha sin cabeza
                'LineWidth',2, ...        % con un ancho de linea 2
                'Marker','.');            % y en el punto (x,y) poner un punto '.'
            
            % la misma flecha girada 180 grados
            quiver(cgx, cgy,...             
                norma.*cos(angulos{i}+pi), norma.*sin(angulos{i}+pi),... 
                esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.');                    
        end            
    end
end
