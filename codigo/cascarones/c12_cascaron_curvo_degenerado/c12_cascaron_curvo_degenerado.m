clear, clc, close all   % borro la memoria, la pantalla y las figuras

% Calculo de los desplazamientos, los angulos de giro, las reacciones, los 
% momentos flectores y las fuerzas cortantes en un cascaron generado

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % algunas constantes que ayudaran en la lectura del codigo
uu = 1; vv = 2; ww = 3; t1 = 4; t2 = 5;

%% Se carga la malla de elementos finitos (LaG, xnod_lo, xnod_up)
global LaG xnod
semiesfera_con_orificio_curvo_S8
%cilindro_scordelli_lo_S8

%% Se calculan los ejes que definen las coordenadas nodales
global vg1 vg2 vg3 t
[vg1, vg2, vg3, t] = vg123_t(xnod_lo, xnod_up);

nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 5*nno;        % numero de grados de libertad globales (5 por nodo)

% nodos vs grados de libertad                     
gdl = reshape(1:ngdl,5,nno)'; % = [(1:5:ngdl)' (2:5:ngdl)' ... (5:5:ngdl)']

%% Se dibuja la malla de elementos finitos
figure
hold on
cg = zeros(nef,3);  % centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y), xnod(LaG(e,[1:8 1]),Z));
   
   % Calculo la posicion del centro de gravedad del elemento finito
   cg(e,:) = mean(xnod(LaG(e,:), :));
   
   h = text(cg(e,X), cg(e,Y), cg(e,Z), num2str(e));
   set(h,'Color', [1 0 0]);
end

% Se indican los nodos restringidos
xyz_apoyo = xnod(restricciones(:,1), :);
plot3(xyz_apoyo(:,X), xyz_apoyo(:,Y), xyz_apoyo(:,Z), 'cx', 'MarkerSize', 40);

% Se indican y numeran el resto de nodos en la malla de EFs
plot3(xnod(:,X), xnod(:,Y), xnod(:,Z), 'b.', 'MarkerSize', 10);
text(xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));

daspect([1 1 1]);
view(3);
grid on;
title('Malla de elementos finitos','FontSize', 26);

%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento rectangular serendipito de 
% 8 nodos junto con sus derivadas con respecto a xi y a eta
% Nforma, dN_dxi, dN_deta
global Nforma
c12_funciones_forma_S8

%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION SELECTIVA)
n_gl_p = 3; % orden de la cuadratura de GL
n_gl_s = 2; % orden de la cuadratura de GL
n_gl_t = 2; % orden de la cuadratura de GL para el espesor (SIEMPRE = 2)

% Calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
[x_gl_p, w_gl_p]  = gausslegendre_quad(n_gl_p);
[x_gl_s, w_gl_s]  = gausslegendre_quad(n_gl_s);
[x_gl_t, w_gl_t]  = gausslegendre_quad(n_gl_t); % integra sobre el espesor

%% matrices constitutivas del elemento en coordenadas locales
Dpp = E/(1-nu^2)* [ 1    nu   0
                    nu   1    0
                    0    0    (1-nu)/2 ];
G = E/(2*(1+nu));    % modulo de rigidez
alpha = 5/6;         % coeficiente de distorsion transversal de la losa de RM
Dsp = diag([alpha*G, alpha*G]);

%% se reserva la memoria RAM de diferentes variables
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);      % grados de libertad de cada elemento finito

% en los siguientes contenedores se almacenara la matriz respectiva para 
% cada punto de integracion: 
Np = cell(nef,n_gl_p,n_gl_p,2); % matrices de funciones de forma calculadas con n_gl_p puntos de integracion
Bp = cell(nef,n_gl_p,n_gl_p,2); % matrices de deformacion generalizada de flexion
Bs = cell(nef,n_gl_s,n_gl_s,2); % matrices de deformacion generalizada de cortante

%% se ensambla la matriz de rigidez global y el vector de fuerzas nodales
%% equivalentes global
for e = 1:nef      % ciclo sobre todos los elementos finitos  
   %% se calcula la matrix de rigidez de flexion Kp del elemento e 
   Kpe = zeros(5*nnoef);
   Mpe = zeros(5*nnoef); % matriz que se utiliza en el calculo de fe
   det_Je_p = zeros(n_gl_p, n_gl_p, 2); % Jacobianos
   
   LaGe = LaG(e,:);     % nodos del EF
   vg1e = vg1(LaGe,:);  % vectores que definen el sistema de coord nodales
   vg2e = vg2(LaGe,:);  % vgi = [nnoef x 3], para i = 1,2,3
   te   = t(LaGe);      % [nnoef x 1] espesor en cada uno de los nodos del EF
  
   for p = 1:n_gl_p
      for q = 1:n_gl_p
         for r = 1:n_gl_t
            xi_gl   = x_gl_p(p);        xe = xnod(LaG(e,:),X);
            eta_gl  = x_gl_p(q);        ye = xnod(LaG(e,:),Y);
            zeta_gl = x_gl_t(r);
            [Bp{e,p,q,r}, det_Je_p(p,q,r)] = Bpp_lam_deg(e, xi_gl, eta_gl, zeta_gl);

            % se arma la matriz de rigidez del elemento e
            Kpe = Kpe + Bp{e,p,q,r}'*Dpp*Bp{e,p,q,r}*det_Je_p(p,q,r)*w_gl_p(p)*w_gl_p(q)*w_gl_t(r);

            % Se evaluan las funciones de forma en los puntos de integracion
            % de Gauss-Legendre
            N = Nforma(xi_gl, eta_gl);

            % Se calcula zb
            zb = zeta_gl*te/2; %(zeta - zeta0)*t/2;               
            
            Np{e,p,q,r} = zeros(3,5*nnoef);
            for i = 1:nnoef
               % Se ensambla la matriz de funciones de forma N
               Ci = [vg1e(i,:)' vg2e(i,:)'];
               Np{e,p,q,r}(:,5*i-4:5*i) = N(i)*[eye(3) -zb(i)*Ci];
            end
                       
            % matriz requerida para calcular el vector de fuerzas nodales 
            % equivalentes (se utiliza la integracion completa)
            Mpe = Mpe + Np{e,p,q,r}'*Np{e,p,q,r}*det_Je_p(p,q,r)*w_gl_p(p)*w_gl_p(q)*w_gl_t(r);
         end
      end
   end     

   
   %% se calcula la matrix de rigidez de flexion Kp del elemento e 
   Kse = zeros(5*nnoef);
   det_Je_s = zeros(n_gl_s, n_gl_s, 2); % Jacobianos   
   for p = 1:n_gl_s
      for q = 1:n_gl_s
         for r = 1:n_gl_t
            xi_gl   = x_gl_s(p);        xe = xnod(LaG(e,:),X);
            eta_gl  = x_gl_s(q);        ye = xnod(LaG(e,:),Y);
            zeta_gl = x_gl_t(r);
            [Bs{e,p,q,r}, det_Je_s(p,q,r)] = Bsp_lam_deg(e, xi_gl, eta_gl, zeta_gl);

            % se arma la matriz de rigidez del elemento e
            Kse = Kse + Bs{e,p,q,r}'*Dsp*Bs{e,p,q,r}*det_Je_s(p,q,r)*w_gl_s(p)*w_gl_s(q)*w_gl_t(r);
         end
      end
   end

   fe = Mpe*tt;
   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:)  gdl(LaG(e,3),:) ...
              gdl(LaG(e,4),:)  gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:)  ];

   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kpe + Kse;
   f(idx{e})        = f(idx{e}) + fe;
end;

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% Se definen las restricciones 
ngdl_res = size(restricciones,1); % numero de grados de libertad restringidos
restric = zeros(ngdl_res,2);
for i = 1:ngdl_res
%                       nodo                direccion           desplazamiento    
   restric(i,:) = [ gdl(restricciones(i,1), restricciones(i,2)) restricciones(i,3) ];
end

%% grados de libertad del desplazamiento conocidos y desconocidos  
c = restric(:,1);   d = setdiff(1:ngdl,c)';

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
ac = restric(:,2);   % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);      % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
aa = nan(ngdl,1);    aa(c) = ac;  aa(d) = ad; % desplazamientos
q  = zeros(ngdl,1);  q(c)  = qd;              % fuerzas nodales equivalentes

%% imprimo los resultados
format short g
disp('Movimientos nodales');
disp('~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(aa,5,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g m, v = %12.4g m, w = %12.4g m, t1 = %12.4g rad, t2 = %12.4g rad\n', ...
           i, vect_mov(i,uu), vect_mov(i,vv), vect_mov(i,ww), vect_mov(i,t1), vect_mov(i,t2));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,5,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0 0 0 0])
      fprintf('Nodo %3d: Fx = %12.4g N, Fy = %12.4g N, Fz = %12.4g N, M1 = %12.4g N-m, M2 = %12.4g N-m\n', ...
         i, q(i,uu), q(i,vv), q(i,ww), q(i,t1), q(i,t2));
   end;
end;

%% Dibujo la malla de elementos finitos y las deformaciones de esta
escala = 2000; % factor de escalamiento de la deformada
xdef   = xnod + escala*vect_mov(:,1:3); % posicion de la deformada

figure; hold on;
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y), xnod(LaG(e,[1:8 1]),Z), 'Color','r');
   line(xdef(LaG(e,[1:8 1]),X), xdef(LaG(e,[1:8 1]),Y), xdef(LaG(e,[1:8 1]),Z), 'Color','b');
end
daspect([1 1 1])
view(3)
grid on
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);
