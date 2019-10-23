clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% ------------------------------------------------------------------------
%% NOTA: este codigo SOLO es apropiado para TENSION PLANA usando elementos
%% rectangulares serendipitos de 8 nodos
%% ------------------------------------------------------------------------

%% DEFINICIÓN DEL PROBLEMA:
% Calcule los desplazamientos y las reacciones en los empotramiento, las
% deformaciones y los esfuerzos de la estructura en TENSION PLANA mostrada 
% en la figura adjunta

%% defino las variables/constantes
X    = 1;           % un par de constantes que ayudaran en la
Y    = 2;           % lectura del codigo
Ee   = 200e9;       % modulo de elasticidad del solido (Pa) = 200 GPa
nue  = 0.30;        % coeficiente de Poisson
te   = 0.01;        % espesor del solido (m)
rhoe = 7850;        % densidad (kg/m^3)
g    = 9.81;        % aceleracion de la gravedad (m/s^2)
be = [0; -rhoe*g];  % vector de fuerzas masicas del elemento

MALLA = 1; % MALLA=1 grafico, MALLA=2 la generada con ANSYS

%% cargar
% xnod - posicion de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
switch MALLA
   case 1
      malla1
   case 2
      malla2
   case 3
      malla3
   otherwise
      error('Malla no especificada')
end
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 2*nno;        % numero de grados de libertad (dos por nodo)
gdl  = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad
nef = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se definen las restricciones 
ngdl_res = size(restricciones,1); % numero de grados de libertad restringidos
restric = zeros(ngdl_res,2);
for i = 1:ngdl_res
%                       nodo                direccion           desplazamiento    
   restric(i,:) = [ gdl(restricciones(i,1), restricciones(i,2)) restricciones(i,3) ];
end

%% Relacion de cargas puntuales
f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
switch MALLA
   case 1
      f(gdl(13,Y)) = -5000;  % carga puntual en el nodo 13 dir Y
      f(gdl(21,Y)) = -5000;  % carga puntual en el nodo 21 dir Y
   case 2
      f(gdl(23,Y)) = -5000;  % carga puntual en el nodo 13 dir Y
      f(gdl(29,Y)) = -5000;  % carga puntual en el nodo 21 dir Y
   case 3  
      % No tenemos cargas puntuales en esta estructura
   otherwise
      error('Malla no especificada')
end  

%% Se definen las cargas distribuidas 
switch MALLA
   case 1
      carga_distr = [];
   case 2
      carga_distr = [];      
   case 3
      %  elem lado  tix   tiy   tjx    tjy  tkx  tky 
      carga_distr = [ ...
         42   123   0         0 0     50000 0    100000  % 42 702 710 717 123
         43   123   0    100000 0    100000 0    100000  % 43 717 727 733 123
         44   123   0    100000 0    100000 0    100000  % 44 733 740 747 123
         45   123   0    100000 0    100000 0    100000  % 45 747 754 760 123
         46   123   0    100000 0    100000 0    100000  % 46 760 767 772 123
         47   123   0    100000 0    100000 0    100000  % 47 772 776 779 123
         48   123   0    100000 0    100000 0    100000  % 48 779 780 778 123
         49   123   0    100000 0    100000 0    100000  % 49 778 774 770 123
         50   123   0    100000 0    100000 0    100000  % 50 770 766 761 123
         51   123   0    100000 0     50000 0       0 ]; % 51 761 757 750 123
   otherwise
      error('Malla no especificada')
end  
nlcd = size(carga_distr,1); % numero de lados con carga distribuida

%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(LaG(e,[1 3 5 7]),X));
   cgy(e) = mean(xnod(LaG(e,[1 3 5 7]),Y));
   h = text(cgx(e), cgy(e), num2str(e)); 
   set(h,'Color', [1 0 0], 'FontSize',16);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'), 'FontSize',16);
axis equal tight
title('Malla de elementos finitos','FontSize',26);

%% Funciones de forma serendipitas del elemento rectangular de 8 nodos:
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa c5_funciones_forma_lagrangianos_rect_2D_8_nodos.m

Nforma = @(xi,eta) [ ...
-((eta - 1)*(xi - 1)*(eta + xi + 1))/4       % N1
((xi^2 - 1)*(eta - 1))/2                     % N2
((eta - 1)*(xi + 1)*(eta - xi + 1))/4        % N3
-((eta^2 - 1)*(xi + 1))/2                    % N4
((eta + 1)*(xi + 1)*(eta + xi - 1))/4        % N5
-((xi^2 - 1)*(eta + 1))/2                    % N6
((eta + 1)*(xi - 1)*(xi - eta + 1))/4        % N7
((eta^2 - 1)*(xi - 1))/2                 ];  % N8

%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
-((eta + 2*xi)*(eta - 1))/4                  % dN1_dxi
eta*xi - xi                                  % dN2_dxi
((eta - 2*xi)*(eta - 1))/4                   % dN3_dxi
1/2 - eta^2/2                                % dN4_dxi
((eta + 2*xi)*(eta + 1))/4                   % dN5_dxi
-xi*(eta + 1)                                % dN6_dxi
-((eta - 2*xi)*(eta + 1))/4                  % dN7_dxi
eta^2/2 - 1/2                            ];  % dN8_dxi

%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
-((2*eta + xi)*(xi - 1))/4                   % dN1_deta
xi^2/2 - 1/2                                 % dN2_deta
((xi + 1)*(2*eta - xi))/4                    % dN3_deta
-eta*(xi + 1)                                % dN4_deta
((2*eta + xi)*(xi + 1))/4                    % dN5_deta
1/2 - xi^2/2                                 % dN6_deta
-((xi - 1)*(2*eta - xi))/4                   % dN7_deta
eta*(xi - 1)                             ];  % dN8_deta

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta
n_gl = 2;                 % orden de la cuadratura de Gauss-Legendre

% El comando:
[x_gl, w_gl]  = gausslegendre_quad(n_gl);
% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
% >> [x_gl,w_gl] = gausslegendre_quad(1)
% x_gl = 0;
% w_gl = 2;
% >> [x_gl,w_gl] = gausslegendre_quad(2)
% x_gl = [  -0.577350269189626;  0.577350269189626 ];
% w_gl = [   1.000000000000000;  1.000000000000000 ];
% >> [x_gl,w_gl] = gausslegendre_quad(3)
% x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
% w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
% >> [x_gl,w_gl] = gausslegendre_quad(4)
% x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
% w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547;
% 0.347854845137453 ];

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl);   % matriz de rigidez global como RALA (sparse)
N = cell(nef,n_gl,n_gl); % contenedor para las matrices de forma
B = cell(nef,n_gl,n_gl); % contenedor para las matrices de deformacion

% matriz constitutiva del elemento para TENSION PLANA
De = [ Ee/(1-nue^2)     Ee*nue/(1-nue^2)  0
       Ee*nue/(1-nue^2) Ee/(1-nue^2)      0
       0                0                 Ee/(2*(1+nue)) ];
idx = cell(nef,1);
for e = 1:nef          % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:) ...
              gdl(LaG(e,3),:)  gdl(LaG(e,4),:) ...
              gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:) ];
   
   % Calculo las matrices de rigidez y el vector de fuerzas nodales
   % equivalentes del elemento
   Ke = zeros(16);
   fe = zeros(16,1);
   det_Je = zeros(n_gl,n_gl); % en esta matriz se almacenaran los Jacobianos

   for p = 1:n_gl
      for q = 1:n_gl
         xi_gl  = x_gl(p);
         eta_gl = x_gl(q);
         
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         NNforma = Nforma(xi_gl, eta_gl);
         
         % Se evaluan las derivadas de las funciones de forma en los puntos
         % de integracion de Gauss-Legendre
         ddN_dxi  = dN_dxi (xi_gl, eta_gl);       xe = xnod(LaG(e,:),X);
         ddN_deta = dN_deta(xi_gl, eta_gl);       ye = xnod(LaG(e,:),Y);
         
         dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
         dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
         
         % Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];
            
         % Se calcula el determinante del Jacobiano
         det_Je(p,q) = det(Je);
         
         N{e,p,q} = zeros(2,2*8);
         B{e,p,q} = zeros(3,2*8);
         for i = 1:8
            % Se ensambla la matriz de funciones de forma N
            N{e,p,q}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                        0           NNforma(i) ];
         
            % Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je(p,q);
            dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je(p,q);
            B{e,p,q}(:,[2*i-1 2*i]) = [ dNi_dx 0          % aqui se ensambla
                                        0      dNi_dy     % y asigna la matriz
                                        dNi_dy dNi_dx ];  % B_i
         end

         % se arma la matriz de rigidez del elemento e
         Ke = Ke + B{e,p,q}'*De*B{e,p,q}*det_Je(p,q)*te*w_gl(p)*w_gl(q);

         % vector de fuerzas nodales equivalentes
         fe = fe + N{e,p,q}'*be*det_Je(p,q)*te*w_gl(p)*w_gl(q);
      end
   end
   
   if any(any(det_Je <= 0))
      error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
   end

   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
end

%% Relacion de las cargas superficiales (vector ft)
ft = sparse(ngdl,1); % fuerzas nodales equivalentes de cargas superficiales
for i = 1:nlcd
   e     = carga_distr(i,1);
   lado  = carga_distr(i,2);
   carga = carga_distr(i,3:8);
   fte = t2ft_R89(xnod(LaG(e,1:8),[X Y]), lado, carga, te);  
   ft(idx{e},:) = ft(idx{e},:) + fte;
end

% Agrego al vector de fuerzas nodales equivalentes las fuerzas
% superficiales calculadas
f = f + ft;

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero', ...
   'FontSize', 26);

%% grados de libertad del desplazamiento conocidos y desconocidos  
c = restric(:,1);   d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = restric(:,2);   % desplazamientos conocidos
qc = zeros(size(d)); % cargas de equilibrio en nodos libres ( = 0 siempre)

%% resuelvo el sistema de ecuaciones
ad = Kdd\((fc+qc)-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd;  q(d) = qc; % fuerzas nodales equivalentes

%% imprimo los resultados
format short g
disp('Nodo   Despl_x (m)   Despl_y (m) = ');     [1:nno; reshape(a,2,nno)]'
disp('Nodo Fuerzas nodales equiv. X, Y (N) = '); [1:nno; reshape(f,2,nno)]'
disp('Nodo Fuerzas nodales equil. X, Y (N) = '); [1:nno; reshape(q,2,nno)]'

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,2,nno)';
escala = 50000;             % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   h1 = line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y)); % original
   set(h1, 'Color', [0 0 1]); % color expresado en notacion RBG entre 0 y 1
   h2 = line(xdef(LaG(e,[1:8 1]),X), xdef(LaG(e,[1:8 1]),Y)); % deformada
   set(h2, 'Color', [1 0 0]);
end
axis equal tight;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside');
title(sprintf('Deformada escalada %d veces',escala), 'FontSize', 26);

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl,n_gl);
esf = cell(nef,n_gl,n_gl);
for e = 1:nef
   ae = a(idx{e});            % desplazamientos de los gdl del elemento e
   
   for pp = 1:n_gl
      for qq = 1:n_gl
         def{e,pp,qq} = B{e,pp,qq}*ae;    % calculo las deformaciones
         esf{e,pp,qq} = De*def{e,pp,qq};  % calculo los esfuerzos
      end
   end
end

%% Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
sx  = zeros(nno,1);
sy  = zeros(nno,1);
sz  = zeros(nno,1);
txy = zeros(nno,1);
txz = zeros(nno,1);
tyz = zeros(nno,1);

ex  = zeros(nno,1);
ey  = zeros(nno,1);
gxy = zeros(nno,1);

A = [ ... 
   3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2
 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4
            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2
 1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4
   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1
 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4
            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2
 3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4 ];

for e = 1:nef
   sx(LaG(e,:),:) = sx(LaG(e,:),:)   + A * [ esf{e,1,1}(1)
                                             esf{e,1,2}(1)
                                             esf{e,2,1}(1)
                                             esf{e,2,2}(1) ];

   sy(LaG(e,:),:) = sy(LaG(e,:),:)   + A * [ esf{e,1,1}(2)
                                             esf{e,1,2}(2)
                                             esf{e,2,1}(2)
                                             esf{e,2,2}(2) ];
                                        
   txy(LaG(e,:),:) = txy(LaG(e,:),:) + A * [ esf{e,1,1}(3)
                                             esf{e,1,2}(3)
                                             esf{e,2,1}(3)
                                             esf{e,2,2}(3) ];                       
                                          
   ex(LaG(e,:),:) = ex(LaG(e,:),:)   + A * [ def{e,1,1}(1)
                                             def{e,1,2}(1)
                                             def{e,2,1}(1)
                                             def{e,2,2}(1) ];

   ey(LaG(e,:),:) = ey(LaG(e,:),:)   + A * [ def{e,1,1}(2)
                                             def{e,1,2}(2)
                                             def{e,2,1}(2)
                                             def{e,2,2}(2) ];
                                        
   gxy(LaG(e,:),:) = gxy(LaG(e,:),:) + A * [ def{e,1,1}(3)
                                             def{e,1,2}(3)
                                             def{e,2,1}(3)
                                             def{e,2,2}(3) ];                                                                 
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;  ex  =  ex./num_elem_ady;
sy  =  sy./num_elem_ady;  ey  =  ey./num_elem_ady;
txy = txy./num_elem_ady;  gxy = gxy./num_elem_ady;

%% Se calculan las deformacion ez en tension plana
ez  = -(nue/Ee)*(sx+sy);

%% Se imprimen y grafican las deformaciones en los nodos
disp('Deformaciones: (Nodo,ex,ey,ez,gxy) = '); 
disp([(1:nno)'  ex  ey  ez  gxy])
figure
subplot(2,2,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ex(LaG(e,:)))
end
ylabel('\epsilon_x','FontSize',26); axis equal tight; colorbar; 

subplot(2,2,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ey(LaG(e,:)))
end
ylabel('\epsilon_y','FontSize',26); axis equal tight; colorbar;

subplot(2,2,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),ez(LaG(e,:)))
end
ylabel('\epsilon_z','FontSize',26); axis equal tight; colorbar;

subplot(2,2,4); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),gxy(LaG(e,:)))
end
ylabel('\gamma_{xy}','FontSize',26); axis equal tight; colorbar;

%% Se imprimen y grafican los esfuerzos en los nodos
disp('Esfuerzos (Pa):  (Nodo,sx,sy,txy) = '); 
disp([(1:nno)'  sx  sy  txy])
figure
subplot(2,2,1); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sx(LaG(e,:)))
end
ylabel('\sigma_x (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(2,2,2); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sy(LaG(e,:)))
end
ylabel('\sigma_y (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(2,2,3); hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),txy(LaG(e,:)))
end
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
disp('Nodo,s1(Pa),s2(Pa),tmax(Pa),angulo(rad) = '); 
disp([(1:nno)'  s1  s2  tmax  ang])
disp('Nodo,Esfuerzos de von Mises (Pa) = ');
disp([(1:nno)'  sv]);

%% s1, s2, taumax
esc = 0.5; % escala para graficar las flechas

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s1(LaG(e,:)))
end

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
title('\sigma_1 (Pa)','FontSize',26); colorbar

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),s2(LaG(e,:)))
end
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
title('\sigma_2 (Pa)','FontSize',26); colorbar

figure;
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),tmax(LaG(e,:)))
end
% Grafique lineas que indiquen direcciones principales de tau_max,
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
title('\tau_{max} (Pa)','FontSize',26); colorbar

figure; hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),sv(LaG(e,:)))
end
ylabel('\sigma_v (Pa)','FontSize',26); axis equal tight; colorbar;
title('Esfuerzos de von Mises (Pa)','FontSize',26);

% Pasando los esfuerzos ya promediados:
export_to_GiD('c5_ejemplo_a',xnod,LaG,a,q,[sx sy sz txy txz tyz]);

% Pasando los puntos de Gauss [RECOMENDADO] !!!
% export_to_GiD('c5_ejemplo_b',xnod,LaG,a,q,esf);                    

%%
return; % bye, bye!
