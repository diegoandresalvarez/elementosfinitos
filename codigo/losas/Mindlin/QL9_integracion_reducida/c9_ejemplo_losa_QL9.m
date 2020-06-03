clear, clc, close all   % borro la memoria, la pantalla y las figuras

% Calculo de los desplazamientos verticales y angulos de giro, las 
% reacciones, los momentos flectores y las fuerzas cortantes en una losa de
% Reissner-Mindlin utilizando los elementos finitos de placa "QL9"

%% defino las variables/constantes
X = 1; Y = 2;        % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
E  = 210e9;          % modulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;            % coeficiente de Poisson
t  = 0.05;           % espesor de la losa (m)
qdistr = -10000;     % carga (N/m^2)

% Definimos la geometria de la losa (generada con generar_malla_losa.m)
load malla_losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nef   = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nnoef = size(LaG,2);  % numero de nodos por EF
nno   = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl  = 3*nno;        % numero de grados de libertad (tres por nodo)

gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del elemento finito
   cgx(e) = mean(xnod(LaG(e,:), X));
   cgy(e) = mean(xnod(LaG(e,:), Y));
   h = text(cgx(e)+0.03, cgy(e)+0.03, num2str(e)); set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Malla de elementos finitos','FontSize', 26);

%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
% junto con sus derivadas con respecto a xi y a eta
c9_funciones_forma_lagrangiano_9_nodos    % Nforma, dN_dxi, dN_deta

%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION SELECTIVA)
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta

% se utilizara integracion SELECTIVA
n_gl_f = 3; % orden de la cuadratura de GL para la integracion de Kf
n_gl_c = 2; % orden de la cuadratura de GL para la integracion de Kc

% El comando:
[x_gl_f, w_gl_f]  = gausslegendre_quad(n_gl_f);
[x_gl_c, w_gl_c]  = gausslegendre_quad(n_gl_c);
% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
% >> [x_gl,w_gl] = gausslegendre_quad(2)
% x_gl = [  -0.577350269189626;  0.577350269189626 ];
% w_gl = [   1.000000000000000;  1.000000000000000 ];
% >> [x_gl,w_gl] = gausslegendre_quad(3)
% x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
% w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];

%% matrices constitutivas del elemento
Df = E/(1-nu^2)* [ 1  nu 0
                   nu 1  0
                   0  0  (1-nu)/2 ];
G = E/(2*(1+nu));  % modulo de rigidez
alpha = 5/6;       % coeficiente de distorsion transversal de la losa de RM
Dc = diag([alpha*G, alpha*G]);
               
Dfg = (t^3/12)*Df; % matriz constitutiva generalizada de flexion
Dcg = t*Dc;        % matriz constitutiva generalizada de cortante

%% se reserva la memoria RAM de diferentes variables
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);      % grados de libertad de cada elemento finito

% en los siguientes contenedores se almacenara la matriz respectiva para 
% cada punto de integracion: 
Nf = cell(nef,n_gl_f,n_gl_f); % matrices de funciones de forma calculadas con n_gl_f puntos de integracion
Bf = cell(nef,n_gl_f,n_gl_f); % matrices de deformacion generalizada de flexion
Bc = cell(nef,n_gl_c,n_gl_c); % matrices de deformacion generalizada de cortante

%% se ensambla la matriz de rigidez global y el vector de fuerzas nodales
%% equivalentes global
for e = 1:nef      % ciclo sobre todos los elementos finitos
   %% se calcula la matrix de rigidez de flexion Kf del elemento e 
   Kfe = zeros(3*nnoef);
   Mfe = zeros(3*nnoef); % matriz que se utiliza en el calculo de fe
   det_Je_f = zeros(n_gl_f); % Jacobianos con n_gl_f puntos de integracion   
   for p = 1:n_gl_f
      for q = 1:n_gl_f
         xi_gl  = x_gl_f(p);        xe = xnod(LaG(e,:),X);
         eta_gl = x_gl_f(q);        ye = xnod(LaG(e,:),Y);
         [Bf{e,p,q}, det_Je_f(p,q)] = Bf_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);

         % se arma la matriz de rigidez del elemento e
         Kfe = Kfe + Bf{e,p,q}'*Dfg*Bf{e,p,q}*det_Je_f(p,q)*w_gl_f(p)*w_gl_f(q);
         
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         Nf{e,p,q} = zeros(3,3*nnoef);
         for i = 1:nnoef
            % Se ensambla la matriz de funciones de forma N
            Nf{e,p,q}(:,3*i-2:3*i) = diag([N(i) N(i) N(i)]);         
         end;                  

         % matriz requerida para calcular el vector de fuerzas nodales 
         % equivalentes (se utiliza la integracion completa)
         Mfe = Mfe + Nf{e,p,q}'*Nf{e,p,q}*det_Je_f(p,q)*w_gl_f(p)*w_gl_f(q);                                              % REVISAR !!!!!!!!!!!!!!!!
      end;
   end; 
       
   %% se calcula la matrix Kc
   Kce = zeros(3*nnoef);   
   det_Je_c = zeros(n_gl_c); % Jacobianos con n_gl_c puntos de integracion
   for p = 1:n_gl_c
      for q = 1:n_gl_c
         xi_gl  = x_gl_c(p);
         eta_gl = x_gl_c(q);         
         xe = xnod(LaG(e,:),X);
         ye = xnod(LaG(e,:),Y);
         [Bc{e,p,q}, det_Je_c(p,q)] = Bc_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);   

         % se arma la matriz de rigidez del elemento e
         Kce = Kce + Bc{e,p,q}'*Dcg*Bc{e,p,q}*det_Je_c(p,q)*w_gl_c(p)*w_gl_c(q);         
      end;
   end;   
   
   %% se calcula el vector de fuerzas nodales equivalentes del elemento e      
   xa = xnod(LaG(e,1),X);   ya = xnod(LaG(e,1),Y);
   xb = xnod(LaG(e,5),X);   yb = xnod(LaG(e,5),Y);
   if (xa >= 0.9999 && xb <= 1.601) && (ya >= 0.9999 && yb <= 2.001)
      ffe = zeros(nnoef, 3); ffe(:,ww) = qdistr;
      ffe = reshape(ffe', 3*nnoef,1);                                                                                             % REVISAR !!!!!!!!!!!!!
   else
      ffe = zeros(3*nnoef,1);
   end;  
   fe = Mfe*ffe;   
   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:)  ...
              gdl(LaG(e,4),:)  gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:)  gdl(LaG(e,9),:) ];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kfe + Kce;
   f(idx{e})        = f(idx{e}) + fe;
end;

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos

% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(abs(xnod(:,X) - 0) < 1e-4);     
lado_y0 = find(abs(xnod(:,Y) - 0) < 1e-4);
lado_x2 = find(abs(xnod(:,X) - 2) < 1e-4);     
lado_y4 = find(abs(xnod(:,Y) - 4) < 1e-4);

c = [ gdl(lado_x0,ww); gdl(lado_x0,ty); 
      gdl(lado_x2,ww); gdl(lado_x2,ty);
      gdl(lado_y0,ww); gdl(lado_y0,tx);
      gdl(lado_y4,ww); gdl(lado_y4,tx) ];

d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = zeros(length(c),1); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = nan(ngdl,1);     a(c) = ac;   a(d) = ad; % desplazamientos
q = zeros(ngdl,1);   q(c) = qd;              % fuerzas nodales equivalentes

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(a,3,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: w = %12.4g m, tx = %12.4g rad, ty = %12.4g rad\n', ...
      i, vect_mov(i,ww), vect_mov(i,tx), vect_mov(i,ty));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,3,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0])
      fprintf('Nodo %3d W = %12.4g N, Mx = %12.4g N-m, My = %12.4g N-m\n', ...
         i, q(i,ww), q(i,tx), q(i,ty));
   end;
end;

%% Se dibuja el plano medio de la malla de elementos finitos y las deformaciones de esta
escala = 5000;            % factor de escalamiento de la deformada
xdef   = escala*vect_mov; % posicion de la deformada
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   fill3(xnod(LaG(e,[1:8 1]),X), ...
         xnod(LaG(e,[1:8 1]),Y), ...
         xdef(LaG(e,[1:8 1]),ww),...
         xdef(LaG(e,[1:8 1]),ww)); %deformada
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
%colorbar('YTick',-0.6:0.05:0)
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
view(3);


%% Se dibuja de la malla de elementos finitos y las deformaciones de esta
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   c9_dibujar_EF_losa(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), ...
      Nforma, a(idx{e}), t, escala, escala);
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
view(3);

%% En los puntos de integracion de Gauss-Legendre calcular:
%% El vector de momentos flectores y torsores (2x2)
%% El vector de fuerzas cortantes (1x1)
sigmag_f = cell(nef, n_gl_f, n_gl_f); % momentos flectores y torsores
sigmag_c = cell(nef, n_gl_c, n_gl_c); % fuerzas cortantes
for e = 1:nef      % ciclo sobre todos los elementos finitos
   for p = 1:n_gl_f
      for q = 1:n_gl_f
         sigmag_f{e,p,q} = Dfg*Bf{e,p,q}*a(idx{e});
      end
   end
   
   for p = 1:n_gl_c
      for q = 1:n_gl_c
         sigmag_c{e,p,q} = Dcg*Bc{e,p,q}*a(idx{e});   
      end
   end
end

%% Se extrapolan los momentos flectores y fuerzas cortantes a los nodos
%% Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);  Qx = zeros(nno,1);
My  = zeros(nno,1);  Qy = zeros(nno,1);
Mxy = zeros(nno,1);

% matriz de extrapolacion de esfuerzos para un elemento lagrangiano de 9
% nodos
A = [ ... 
   3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2
 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4
            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2
 1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4
   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1
 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4
            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2
 3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4
             1/4,             1/4,             1/4,             1/4 ];

for e = 1:nef
   Mx(LaG(e,:),:)  = Mx(LaG(e,:),:)  + A * [ sigmag_f{e,1,1}(1)
                                             sigmag_f{e,1,2}(1)
                                             sigmag_f{e,2,1}(1)
                                             sigmag_f{e,2,2}(1) ];

   My(LaG(e,:),:)  = My(LaG(e,:),:)  + A * [ sigmag_f{e,1,1}(2)
                                             sigmag_f{e,1,2}(2)
                                             sigmag_f{e,2,1}(2)
                                             sigmag_f{e,2,2}(2) ];
                                        
   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ sigmag_f{e,1,1}(3)
                                             sigmag_f{e,1,2}(3)
                                             sigmag_f{e,2,1}(3)
                                             sigmag_f{e,2,2}(3) ];
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
Mx  =  Mx./num_elem_ady;
My  =  My./num_elem_ady;
Mxy = Mxy./num_elem_ady;

%% Se imprimen y grafican los esfuerzos en los nodos
disp('Esfuerzos (Pa):  (Nodo,Mx,My,Mxy) = '); 
disp([(1:nno)'  Mx  My  Mxy])
figure
subplot(1,3,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),Mx(LaG(e,1:8)))
end;
ylabel('Momentos Mx (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')

subplot(1,3,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),My(LaG(e,1:8)))
end;
ylabel('Momentos My (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')

subplot(1,3,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),Mxy(LaG(e,1:8)))
end;
ylabel('Momentos Mxy (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')


sx = Mx; sy = My; txy = Mxy;
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

%% imprimo los resultados
disp('Nodo, M1(N-m/m), M2(N-m/m),Mtorsor_max(N-m/m), angulo(rad) = '); 
disp([(1:nno)'  s1  s2  tmax  ang])

%% s1, s2, taumax
esc = 0.5; % escala para graficar las flechas

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),s1(LaG(e,1:8)))
end;

% Grafique lineas que indican las direcciones principales de sigma_1
norma = 1; % = s1 si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...   % En el nodo grafique una flecha (linea)
   norma.*cos(ang),norma.*sin(ang),... % indicando la direccion principal de sigma_1
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(xnod(:,X),xnod(:,Y),...   % la misma flecha ahora en la otra direccion,
   norma.*cos(ang+pi),norma.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('M_1 (N-m/m)','FontSize',26); colorbar

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),s2(LaG(e,1:8)))
end;
% Grafique lineas que indiquen direcciones principales de sigma_2
norma = 1; % = s2 si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...             % flecha indicando la direccion
   norma.*cos(ang+pi/2),norma.*sin(ang+pi/2),... % principal de sigma_2
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
   norma.*cos(ang-pi/2),norma.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('M_2 (N-m/m)','FontSize',26); colorbar

figure;
hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),tmax(LaG(e,1:8)))
end;
% Grafique lineas que indiquen direcciones principales de Mtorsor_max,
norma = 1; % = tmax si quiere proporcional
quiver(xnod(:,X),xnod(:,Y), ...
       norma.*cos(ang+pi/4),norma.*sin(ang+pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang-pi/4),norma.*sin(ang-pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang+3*pi/4),norma.*sin(ang+3*pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang-3*pi/4),norma.*sin(ang-3*pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('Mtorsor_{max} (N-m/m)','FontSize',26); colorbar

%% FALTA DIAGRAMA DE LAS FUERZAS CORTANTES


%% Refuerzo óptimo de la losa
esc = 0.5; % escala para graficar las flechas

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),s1(LaG(e,1:8)))
end;

% Grafique lineas que indican las direcciones principales de sigma_1
norma = double(s1>0);
quiver(xnod(:,X),xnod(:,Y),...   % En el nodo grafique una flecha (linea)
   norma.*cos(ang),norma.*sin(ang),... % indicando la direccion principal de sigma_1
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(xnod(:,X),xnod(:,Y),...   % la misma flecha ahora en la otra direccion,
   norma.*cos(ang+pi),norma.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('M_1 (N-m/m): Refuerzo optimo parte superior de la losa','FontSize',26); colorbar

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,1:8),X),xnod(LaG(e,1:8),Y),s2(LaG(e,1:8)))
end;
% Grafique lineas que indiquen direcciones principales de sigma_2
norma = double(s2<0);
quiver(xnod(:,X),xnod(:,Y),...             % flecha indicando la direccion
   norma.*cos(ang+pi/2),norma.*sin(ang+pi/2),... % principal de sigma_2
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
   norma.*cos(ang-pi/2),norma.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('M_2 (N-m/m): Refuerzo optimo parte inferior de la losa','FontSize',26); colorbar

return; % bye, bye!
