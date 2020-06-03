%% Programa que ilustra el modo de energia nula que aparece en el EF QL9
%  cuando se hace integracion selectiva

% Calculo de los desplazamientos verticales y angulos de giro, las 
% reacciones, los momentos flectores y las fuerzas cortantes en una losa de
% Mindlin utilizando los elementos finitos de placa "QL9"

%%
clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
X = 1; Y = 2;        % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
E  = 210e9;          % modulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;            % coeficiente de Poisson
t  = 0.05;           % espesor de la losa (m)
qdistr = -10000;     % carga (N/m^2)

% Definimos la geometria de la losa (creada con "generar_malla_losa.m")
load malla_losa_MEN
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nef   = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nnoef = size(LaG,2);  % numero de nodos por EF
nno   = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl  = 3*nno;        % numero de grados de libertad (tres por nodo)

gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%% Se dibuja la malla de elementos finitos
figure; 
hold on;
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del elemento finito
   cgx = mean(xnod(LaG(e,:), X));
   cgy = mean(xnod(LaG(e,:), Y));
   text(cgx+0.03, cgy+0.03, num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'rx');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis([-0.5, 2.5, -0.5, 4.5])
title('Malla de una losa con EFs QL9');

%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
% junto con sus derivadas con respecto a xi y a eta
funciones_forma_lagrangiano_9_nodos    % Nforma, dN_dxi, dN_deta

%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION SELECTIVA)
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta

% se utilizara integracion COMPLETA
%{
n_gl_b = 3; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 3; % orden de la cuadratura de GL para la integracion de Ks
%}

% se utilizara integracion SELECTIVA
n_gl_b = 3; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 2; % orden de la cuadratura de GL para la integracion de Ks

% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
[x_gl_b, w_gl_b]  = gausslegendre_quad(n_gl_b);
[x_gl_s, w_gl_s]  = gausslegendre_quad(n_gl_s);

%% matrices constitutivas del elemento
Db = (E/(1-nu^2))* [ 1    nu   0
                     nu   1    0
                     0    0    (1-nu)/2 ];              
G = E/(2*(1+nu));  % modulo de rigidez
alpha = 5/6;       % coeficiente de distorsion transversal de la losa de RM
Ds = diag([alpha*G, alpha*G]);
               
Dbg = (t^3/12)*Db; % matriz constitutiva generalizada de flexion
Dsg = t*Ds;        % matriz constitutiva generalizada de cortante

%% se reserva la memoria RAM de diferentes variables
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);      % grados de libertad de cada elemento finito

% en los siguientes contenedores se almacenara la matriz respectiva para 
% cada punto de integracion: 
NN = cell(nef,n_gl_b,n_gl_b); % matrices de funciones de forma calculadas con n_gl_b puntos de integracion
Bb = cell(nef,n_gl_b,n_gl_b); % matrices de deformacion generalizada de flexion
Bs = cell(nef,n_gl_s,n_gl_s); % matrices de deformacion generalizada de cortante

%% se ensambla la matriz de rigidez global K y el vector de fuerzas nodales
%% equivalentes global f
for e = 1:nef      % ciclo sobre todos los elementos finitos
   xe = xnod(LaG(e,:),X);   
   ye = xnod(LaG(e,:),Y);    
    
   %% se calcula la matrix de rigidez de flexion Kb del elemento e 
   Kbe = zeros(3*nnoef);
   % Mbe = zeros(3*nnoef); % matriz que se utiliza en el calculo de fe   
   det_Je_b = zeros(n_gl_b); % Jacobianos con n_gl_b puntos de integracion   
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = x_gl_b(q);
         [Bb{e,p,q}, det_Je_b(p,q)] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);

         % se arma la matriz de rigidez del elemento e
         Kbe = Kbe + Bb{e,p,q}'*Dbg*Bb{e,p,q}*det_Je_b(p,q)*w_gl_b(p)*w_gl_b(q);
      end
   end
   
   %% se calcula la matrix Ks
   Kse = zeros(3*nnoef);   
   det_Je_s = zeros(n_gl_s); % Jacobianos con n_gl_s puntos de integracion
   for p = 1:n_gl_s
      for q = 1:n_gl_s
         xi_gl  = x_gl_s(p);        
         eta_gl = x_gl_s(q);
         [Bs{e,p,q}, det_Je_s(p,q)] = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);   

         % se arma la matriz de rigidez del elemento e
         Kse = Kse + Bs{e,p,q}'*Dsg*Bs{e,p,q}*det_Je_s(p,q)*w_gl_s(p)*w_gl_s(q);         
      end
   end 
   
   %{
   %% se calcula la matriz NN
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b(p);
         eta_gl = x_gl_b(q);
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         % Se ensambla la matriz de funciones de forma N
         NN{e,p,q} = zeros(3,3*nnoef);
         for i = 1:nnoef            
            NN{e,p,q}(:,3*i-2:3*i) = diag([N(i) N(i) N(i)]);
         end
   
         % matriz requerida para calcular el vector de fuerzas nodales 
         % equivalentes (se utiliza la integracion completa)
         Mbe = Mbe + NN{e,p,q}'*NN{e,p,q}*det_Je_b(p,q)*w_gl_b(p)*w_gl_b(q);                                              % REVISAR !!!!!!!!!!!!!!!!   
      end
   end  
   %}
   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:)  ...
              gdl(LaG(e,4),:)  gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:)  gdl(LaG(e,9),:) ];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kbe + Kse;
   % f(idx{e})        = f(idx{e}) + fe;
end

f(gdl(45,ww)) = -10;

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos
c = [ gdl( 1,ww) gdl( 1,tx) gdl( 1,ty) ...
      gdl(37,ww) gdl(37,tx) gdl(37,ty) ];

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
vect_mov = reshape(a,3,nno)'; % vector de movimientos
%{
format short g
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nno
   fprintf('Nodo %3d: w = %12.4g m, tx = %12.4g rad, ty = %12.4g rad\n', ...
      i, vect_mov(i,ww), vect_mov(i,tx), vect_mov(i,ty));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,3,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0])
      fprintf('Nodo %3d W = %12.4g N, Mx = %12.4g N-m, My = %12.4g N-m\n', ...
         i, q(i,ww), q(i,tx), q(i,ty));
   end
end
%}

%% Se dibuja el plano medio de la malla de elementos finitos y las deformaciones de esta
escala = 5000;            % factor de escalamiento de la deformada
%{
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
%}

%% Se dibuja de la malla de elementos finitos y las deformaciones de esta
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   dibujar_EF_QL9(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), ...
      Nforma, a(idx{e}), t, escala, escala);
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
view(3);

%% En los puntos de integracion de Gauss-Legendre calcular:
%% El vector de momentos flectores y torsores (2x2)
%% El vector de fuerzas cortantes (1x1)
[x_gl_b, w_gl_b]  = gausslegendre_quad(2);
[x_gl_s, w_gl_s]  = gausslegendre_quad(1);

%% se calcula de nuevo Bb y Bs en cada punto de GL
Bb = cell(nef,2,2); % matrices de deformacion generalizada de flexion
Bs = cell(nef);     % matrices de deformacion generalizada de cortante
for e = 1:nef      % ciclo sobre todos los elementos finitos
    xe = xnod(LaG(e,:),X);
    ye = xnod(LaG(e,:),Y);    
    
    %% se calcula la matrix Bb en los puntos de integracion de GL para el 
    % calculo de los momentos flectores y torsores
    for p = 1:2
      for q = 1:2
         xi_gl  = x_gl_b(p);
         eta_gl = x_gl_b(q);
         Bb{e,p,q} = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta);
      end
   end
       
   %% se calcula la matrix Bs en los puntos de integracion de GL para el 
   % calculo de las fuerzas cortantes
   xi_gl  = x_gl_s(1);
   eta_gl = x_gl_s(1);
   Bs{e} = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta);
end

%% Se calculan los momentos y las fuerzas en los puntos de GL
sigmag_b = cell(nef, 2, 2); % momentos flectores y torsores
sigmag_s = cell(nef, 1);    % fuerzas cortantes
for e = 1:nef               % ciclo sobre todos los elementos finitos
   for p = 1:2
      for q = 1:2
         sigmag_b{e,p,q} = Dbg*Bb{e,p,q}*a(idx{e});
      end
   end
   
   sigmag_s{e} = Dsg*Bs{e}*a(idx{e});
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
   Mx(LaG(e,:),:)  = Mx(LaG(e,:),:)  + A * [ sigmag_b{e,1,1}(1)
                                             sigmag_b{e,1,2}(1)
                                             sigmag_b{e,2,1}(1)
                                             sigmag_b{e,2,2}(1) ];

   My(LaG(e,:),:)  = My(LaG(e,:),:)  + A * [ sigmag_b{e,1,1}(2)
                                             sigmag_b{e,1,2}(2)
                                             sigmag_b{e,2,1}(2)
                                             sigmag_b{e,2,2}(2) ];

   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ sigmag_b{e,1,1}(3)
                                             sigmag_b{e,1,2}(3)
                                             sigmag_b{e,2,1}(3)
                                             sigmag_b{e,2,2}(3) ];

   Qx(LaG(e,:),:) = Qx(LaG(e,:),:) + sigmag_s{e}(1);
   Qy(LaG(e,:),:) = Qy(LaG(e,:),:) + sigmag_s{e}(2);

   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
Mx  =  Mx./num_elem_ady;
My  =  My./num_elem_ady;
Mxy = Mxy./num_elem_ady;
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady;

%% Se grafican los momentos
figure
subplot(1,3,1); plot_M_or_Q(nef, xnod, LaG, Mx,  'Momentos Mx (N-m/m)');
subplot(1,3,2); plot_M_or_Q(nef, xnod, LaG, My,  'Momentos My (N-m/m)');
subplot(1,3,3); plot_M_or_Q(nef, xnod, LaG, Mxy, 'Momentos Mxy (N-m/m)');

%% Se grafican los cortantes
figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Qx,  'Cortantes Qx (N/m)');
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Qy,  'Cortantes Qy (N/m)');

%% Se calculan y grafican para cada elemento los momentos principales y
%% sus direcciones
Mt_max = sqrt(((Mx-My)/2).^2 + Mxy.^2); % momento torsion maximo
Mf1_xy = (Mx+My)/2 + Mt_max;            % momento flector maximo
Mf2_xy = (Mx+My)/2 - Mt_max;            % momento flector minimo
ang  = 0.5*atan2(2*Mxy, Mx-My);         % angulo de inclinacion de Mf1_xy

%% Mf1_xy, Mf2_xy, Mt_max
figure
subplot(1,3,1); plot_M_or_Q(nef, xnod, LaG, Mf1_xy, 'Mf1_{xy} (N-m/m)', { ang })
subplot(1,3,2); plot_M_or_Q(nef, xnod, LaG, Mf2_xy, 'Mf2_{xy} (N-m/m)', { ang+pi/2 })
subplot(1,3,3); plot_M_or_Q(nef, xnod, LaG, Mt_max, 'Mt_{max} (N-m/m)', { ang+pi/4, ang-pi/4 })

%% Se calculan y grafican los cortantes maximos, junto con su angulo de inclinacion
Q_max = hypot(Qx, Qy);
ang   = atan2(Qy, Qx);

figure
plot_M_or_Q(nef, xnod, LaG, Q_max, 'Q_{max} (N/m)', { ang })

%%
return; % bye, bye!

%%
function plot_M_or_Q(nef, xnod, LaG, variable, texto, angulos)
    X = 1; Y = 2;
    hold on; 
    colorbar;
    % Por simplicidad no se graficaran los resultados asociados al nodo 9
    for e = 1:nef  
       fill(xnod(LaG(e,1:8),X), xnod(LaG(e,1:8),Y), variable(LaG(e,1:8)));
    end
    axis equal tight
    colormap jet
    title(texto, 'FontSize',20);
   
    esc = 0.5;
    if nargin == 6
        norma = 1; % = variable % si se quiere proporcional
        for i = 1:length(angulos)
            % se indica la flecha de la direccion principal
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}), norma.*sin(angulos{i}),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off',... % una flecha sin cabeza
                'LineWidth',2, ...        % con un ancho de linea 2
                'Marker','.');            % y en el punto (x,y) poner un punto '.'
            
            % la misma flecha girada 180 grados
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}+pi), norma.*sin(angulos{i}+pi),... 
                esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.');                    
        end            
    end
end