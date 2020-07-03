clear, clc, close all   % borro la memoria, la pantalla y las figuras

% Calculo de los desplazamientos, los angulos de giro, las reacciones, los 
% momentos flectores y las fuerzas cortantes en un cascaron de 
% Reissner-Mindlin utilizando los elementos finitos QLLL

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % algunas constantes que ayudaran en la lectura del codigo
uu = 1; vv = 2; ww = 3; tx = 4; ty = 5; tz = 6;

%% Se carga la malla de elementos finitos
semiesfera_con_orificio_Q4
%cilindro_Q4

nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 6*nno;        % numero de grados de libertad globales (seis por nodo)

% nodos vs grados de libertad                     
gdl = reshape(1:ngdl,6,nno)'; % = [(1:6:ngdl)' (2:6:ngdl)' ... (6:6:ngdl)']

%% Se dibuja la malla de elementos finitos
figure; hold on;
cg = zeros(3,nef);  % centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:4 1]),X), xnod(LaG(e,[1:4 1]),Y), xnod(LaG(e,[1:4 1]),Z));
   
   % Calculo la posicion del centro de gravedad del elemento finito
   cg(:,e) = mean(xnod(LaG(e,:), :));
   
   h = text(cg(X,e), cg(Y,e), cg(Z,e), num2str(e), 'Color', [1 0 0]);
end
plot3(xnod(:,X), xnod(:,Y), xnod(:,Z), 'r*');
text (xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));
daspect([1 1 1]);
view(3);
grid on;
title('Malla de elementos finitos', 'FontSize', 26);

%% Se cargan las funciones de forma junto con sus derivadas
% Se cargan las funciones de forma del elemento rectangular de 4 nodos 
% junto con sus derivadas con respecto a xi y a eta
% Nforma, dN_dxi, dN_deta
funciones_forma_Q4

%% parametros de la cuadratura de Gauss-Legendre (INTEGRACION COMPLETA)
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta

% se utilizara integracion COMPLETA (n_gl = n_gl_f = n_gl_c = 2)
n_gl = 2; % orden de la cuadratura de GL

% Calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% matrices constitutivas del elemento en coordenadas locales
Dpp = E/(1-nu^2) * [ 1  nu 0
                     nu 1  0
                     0  0  (1-nu)/2 ];
G = E/(2*(1+nu));    % modulo de rigidez
alpha = 5/6;         % coeficiente de distorsion transversal de la losa de RM
Dsp = diag([alpha*G, alpha*G]);

Dmgp = t*Dpp;        % matriz constitutiva generalizada de membrana
Dbgp = (t^3/12)*Dpp; % matriz constitutiva generalizada de flexion
Dsgp = t*Dsp;        % matriz constitutiva generalizada de cortante

%% se reserva la memoria RAM de diferentes variables
K   = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);      % grados de libertad de cada EF
T   = cell(nef, 1);      % matriz de transformacion de cada EF
lambda = cell(nef, 1);   % matriz de transformacion de cada EF

% en los siguientes contenedores se almacenara la matriz respectiva para 
% cada punto de integracion: 
NN = cell(nef,n_gl,n_gl); % matrices de funciones de forma calculadas con n_gl puntos de integracion
Bb = cell(nef,n_gl,n_gl); % matrices de deformacion generalizada de flexion
Bm = cell(nef,n_gl,n_gl); % matrices de deformacion generalizada de membrana
Bs = cell(nef,n_gl,n_gl); % matrices de deformacion generalizada de cortante

%% se ensambla la matriz de rigidez global y el vector de fuerzas nodales
%% equivalentes global
for e = 1:nef      % ciclo sobre todos los elementos finitos   
   %% se calcula la matriz de transformacion de coordenadas (sec 10.7.1)
   % ALGORITMO 1:
     [T{e}, lambda{e}] = calculo_T(xnod(LaG(e,[1 2 3]),:));
   %% se calcula la matriz de transformacion de coordenadas (sec 10.7.2)
   % ALGORITMO 2:
   % [T{e}, lambda{e}] = calculo_T2(xnod(LaG(e,[1 2 3]),:));   
     
   %% se convierten las coordenadas globales a locales
   xnod_e_loc = xnod(LaG(e,:),:)*lambda{e}'; % eq 8.51 -> xpT = xT*lambdaeT
   xe = xnod_e_loc(:,X);
   ye = xnod_e_loc(:,Y);
   ze = xnod_e_loc(:,Z);
   
   % se verifica que todos los puntos esten sobre el plano x'y'
   if range(ze) > 1e-3
      error('Los nodos del elemento %d no son coplanares\n', e);
   end
   
   %% se separa la memoria RAM
   Ke     = zeros(6*nnoef);
   Mloc   = zeros(5*nnoef);  % matriz que se utiliza en el calculo de fe
   det_Je = zeros(n_gl);     % Jacobianos con n_gl puntos de integracion
   
   %% se calcula Ke
   for p = 1:n_gl
      for q = 1:n_gl
         xi_gl  = x_gl(p);
         eta_gl = x_gl(q);
         
         [Bb{e,p,q}, Bm{e,p,q}, Bs{e,p,q}, det_Je(p,q)] = ...
                                Bb_Bm_Bs_QLLL(xi_gl, eta_gl, xe, ye, T{e});
         
         Kbe = Bb{e,p,q}'*Dbgp*Bb{e,p,q}*det_Je(p,q)*w_gl(p)*w_gl(q);
         Kme = Bm{e,p,q}'*Dmgp*Bm{e,p,q}*det_Je(p,q)*w_gl(p)*w_gl(q);        
         Kse = Bs{e,p,q}'*Dsgp*Bs{e,p,q}*det_Je(p,q)*w_gl(p)*w_gl(q);         

         Ke = Ke + (Kbe + Kse + Kme);
        
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         NN{e,p,q} = zeros(5,5*nnoef);
         for i = 1:nnoef
            % Se ensambla la matriz de funciones de forma N
            NN{e,p,q}(:,5*i-4:5*i) = diag([N(i) N(i) N(i) N(i) N(i)]);
         end

         % matriz requerida para calcular el vector de fuerzas nodales 
         % equivalentes (se utiliza la integracion completa)
         Mloc = Mloc + NN{e,p,q}'*NN{e,p,q}*det_Je(p,q)*w_gl(p)*w_gl(q);
      end
   end     

   %% se calcula el vector de fuerzas nodales equivalentes del elemento e
   fe = T{e}'*Mloc*T{e} * tt;
   
   %% se asocian los grados de libertad del elemento locales a los globales
   idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:) ];
   
   %% se procede al ensamblaje
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:) + fe;
end

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
disp('Movimientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(aa,6,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g m, v = %12.4g m, w = %12.4g m, tx = %12.4g rad, ty = %12.4g rad, tz = %12.4g rad\n', ...
           i, vect_mov(i,uu), vect_mov(i,vv), vect_mov(i,ww), vect_mov(i,tx), vect_mov(i,ty), vect_mov(i,tz));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,6,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0 0 0 0])
      fprintf('Nodo %3d: Fx = %12.4g N, Fy = %12.4g N, Fz = %12.4g N, Mx = %12.4g N-m, My = %12.4g N-m, Mz = %12.4g N-m\n', ...
         i, q(i,uu), q(i,vv), q(i,ww), q(i,tx), q(i,ty), q(i,tz));
   end
end

%% Dibujo la malla de elementos finitos y las deformaciones de esta
escala = 100; % factor de escalamiento de la deformada
xdef   = xnod + escala*vect_mov(:,1:3); % posicion de la deformada

figure; hold on;
for e = 1:nef
   line(xnod(LaG(e,[1:4 1]),X), xnod(LaG(e,[1:4 1]),Y), xnod(LaG(e,[1:4 1]),Z), 'Color','r');
   line(xdef(LaG(e,[1:4 1]),X), xdef(LaG(e,[1:4 1]),Y), xdef(LaG(e,[1:4 1]),Z), 'Color','b');
end
daspect([1 1 1]);
view(3);
grid on;
title(sprintf('Deformada escalada %d veces', escala), 'FontSize', 20);



% FALTA HACER LO QUE SIGUE:
% Hacer la conversion de los esfuerzos de ejes locales a globales mediante
% la ecuacion "sigma = T'*sigmap*T"  
%
% Tenga en cuenta que no tiene sentido graficar los momentos en coordenadas
% locales, ya que la eleccion de los ejes locales es algo arbitraria y 
% depende de la numeracion de la malla. Por lo tanto solo se deben reportar
% los graficos de los momentos maximos/minimos

% El diagrama de momento solo no sirve... hay que incluirle las fuerzas de
% membrana y calcular en la superficie superior e inferior los esfuerzos
% minimos y maximos

return

%% En los puntos de integracion de Gauss-Legendre calcular en coord locales: 
%% El vector de momentos flectores y torsores (1x1)
n_gl_f = n_gl;
%% El vector de fuerzas cortantes (1x1)
n_gl_c = n_gl;

sigmag_f = cell(nef, n_gl_f, n_gl_f); % momentos flectores y torsores
sigmag_c = cell(nef, n_gl_c, n_gl_c); % fuerzas cortantes

mom_f = cell(nef, n_gl_f, n_gl_f); % momentos flectores y torsores minimos/maximos

for e = 1:nef       % ciclo sobre todos los elementos finitos
   TT = lambda{e}'; % Tenga en cuenta que lambda = T' (ver main.pdf)  
   for p = 1:n_gl_f
      for q = 1:n_gl_f
         sigmagp_f_e_p_q = Dbgp*Bb{e,p,q}*aa(idx{e});
         
         mom_f{e,p,q} = TT*[sigmagp_f_e_p_q(1) sigmagp_f_e_p_q(3) 0
                            sigmagp_f_e_p_q(3) sigmagp_f_e_p_q(2) 0
                            0                   0                 0 ]*TT';
                        
         [evec, eval] = eig(mom_f{e,p,q})
      end
   end
   
   for p = 1:n_gl_c
      for q = 1:n_gl_c
         sigmagp_m{e,p,q} = Dmgp*Bm{e,p,q}*aa(idx{e});
         sigmagp_c{e,p,q} = Dsgp*Bs{e,p,q}*aa(idx{e});             
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
   3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2   % 1
            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2   % 2
   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1   % 3
            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2 ];% 4

for e = 1:nef
   Mx(LaG(e,:),:)  = Mx(LaG(e,:),:)  + A * [ sigmagp_f{e,1,1}(1)
                                             sigmagp_f{e,1,2}(1)
                                             sigmagp_f{e,2,1}(1)
                                             sigmagp_f{e,2,2}(1) ];

   My(LaG(e,:),:)  = My(LaG(e,:),:)  + A * [ sigmagp_f{e,1,1}(2)
                                             sigmagp_f{e,1,2}(2)
                                             sigmagp_f{e,2,1}(2)
                                             sigmagp_f{e,2,2}(2) ];
                                        
   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ sigmagp_f{e,1,1}(3)
                                             sigmagp_f{e,1,2}(3)
                                             sigmagp_f{e,2,1}(3)
                                             sigmagp_f{e,2,2}(3) ];
                                          
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
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),Mx(LaG(e,1:4)))
end;
ylabel('Momentos Mx (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')
view(3); grid on

subplot(1,3,2); hold on;
for e = 1:nef
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),My(LaG(e,1:4)))
end;
ylabel('Momentos My (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')
view(3); grid on

subplot(1,3,3); hold on;
for e = 1:nef
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),Mxy(LaG(e,1:4)))
end;
ylabel('Momentos Mxy (N-m/m)','FontSize',26); axis equal tight;
colorbar('Location','SouthOutside')
view(3); grid on

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
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),s1(LaG(e,1:4)))
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
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),s2(LaG(e,1:4)))
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
   fill3(xnod(LaG(e,1:4),X),xnod(LaG(e,1:4),Y),xnod(LaG(e,1:4),Z),tmax(LaG(e,1:4)))
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

return; % bye, bye!

%% FALTA DIAGRAMA DE LAS FUERZAS CORTANTES
