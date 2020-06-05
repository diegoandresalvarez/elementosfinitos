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

% se utilizara integracion REDUCIDA (NO FUNCIONA)
%{
n_gl_b = 2; % orden de la cuadratura de GL para la integracion de Kb
n_gl_s = 2; % orden de la cuadratura de GL para la integracion de Ks
%}

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

   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:)  ...
              gdl(LaG(e,4),:)  gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:)  gdl(LaG(e,9),:) ];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kbe + Kse;
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

%% Se dibuja el plano medio de la malla de elementos finitos y las deformaciones de esta
escala = 5000;            % factor de escalamiento de la deformada
%{
xdef     = escala*vect_mov; % posicion de la deformada
vect_mov = reshape(a,3,nno)'; % vector de movimientos

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
colormap jet
view(3);
%}

%% Se dibuja de la malla de elementos finitos y las deformaciones de esta
figure; 
hold on; 
grid on;
colorbar
for e = 1:nef
   dibujar_EF_Q89_RM(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), ...
      Nforma, a(idx{e}), t, escala, escala);
end
daspect([1 1 1]); % similar a "axis equal", pero en 3D
axis tight
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
colormap jet
view(3);

%% bye, bye !!!
return