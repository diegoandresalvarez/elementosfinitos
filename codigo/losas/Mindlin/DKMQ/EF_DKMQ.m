% Calculo de los desplazamientos en una placa utilizando la teoria de
% Reissner-Mindlin y el elemento finito de placa DKMQ 
% Por:
% Diego Andres Alvarez Marin
% Sebastian Jaramillo Moreno

clear, clc, close all % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

E  = 210e9;       % modulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;         % coeficiente de Poisson
t  = 0.05;        % espesor de la losa (m)
q  = -10000;      % carga (N/m^2)

% Definimos la geometria de la losa
losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1 2 3 4 1]),X), xnod(LaG(e,[1 2 3 4 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X))/2;
   cgy(e) = (xnod(LaG(e,2),Y) + xnod(LaG(e,3),Y))/2;
   text(cgx(e), cgy(e), num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
title('Malla de elementos finitos');

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta
n_gl = 2;                 % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% Se leen las funciones de forma
funciones_de_forma;

%% matrices constitutivas
De = (E/(1-nu^2)) * [ 1  nu 0
                      nu 1  0
                      0  0  (1-nu)/2 ];
               
Dbe = (t^3/12)*De;         % matriz constitutiva de flexion generalizada  

G = E/(2*(1+nu));          % modulo de cortante
Dse = (5/6)*G*t*eye(2);    % matriz constitutiva de cortante generalizada                       

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
idx = cell(nef, 1);    % grados de libertad de cada elemento finito

for e = 1:nef      % ciclo sobre todos los elementos finitos
    idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:) ];
               
    % Calculo las matrices de rigidez y el vector de fuerzas nodales
    % equivalentes del elemento
    Kbe = zeros(12);
    Kse = zeros(12);
    fe  = zeros(12,1);
    det_JeN = zeros(n_gl,n_gl); % en esta matriz se almacenaran los Jacobianos
    det_JeP = zeros(n_gl,n_gl);
    
    for pp = 1:n_gl
        for qq = 1:n_gl
            xi_gl  = x_gl(pp);
            eta_gl = x_gl(qq);

            % Se evaluan las funciones de forma en los puntos de integracion
            % de Gauss-Legendre
            NNforma = Nforma(xi_gl, eta_gl);
            PPforma = Pforma(xi_gl, eta_gl);
            
            % Se evaluan las derivadas de las funciones de forma en los puntos
            % de integracion de Gauss-Legendre
            ddN_dxi  = dN_dxi (xi_gl, eta_gl);       xe = xnod(LaG(e,:),X);
            ddN_deta = dN_deta(xi_gl, eta_gl);       ye = xnod(LaG(e,:),Y);
            ddP_dxi  = dP_dxi (xi_gl, eta_gl);       
            ddP_deta = dP_deta(xi_gl, eta_gl);
            
            xee = [xe; xe(1)];    yee = [ye; ye(1)];
            
            % Coordenadas de los nodos secundarios
            xes = xe + diff(xee)/2;         yes = ye + diff(yee)/2;
            
            dx_dxi_N  = sum(ddN_dxi  .* xe);    dy_dxi_N  = sum(ddN_dxi  .* ye);
            dx_deta_N = sum(ddN_deta .* xe);    dy_deta_N = sum(ddN_deta .* ye);
            dx_dxi_P  = sum(ddP_dxi  .* xe);    dy_dxi_P  = sum(ddP_dxi  .* ye);
            dx_deta_P = sum(ddP_deta .* xe);    dy_deta_P = sum(ddP_deta .* ye);
            
            % Se ensambla la matriz Jacobiana del elemento
            % Matriz Jacobiana de funciones principales
            JeN = [ dx_dxi_N    dy_dxi_N
                    dx_deta_N   dy_deta_N ];
            JeP = [ dx_dxi_P    dy_dxi_P
                    dx_deta_P   dy_deta_P ];
                    
            % Se calcula el determinante del Jacobiano
            det_JeN(pp,qq) = det(JeN);
            det_JeP(pp,qq) = det(JeP);
            
            N{e,pp,qq}     = zeros(3,12);
            Bb_t{e,pp,qq}  = zeros(3,12);
            Bb_dt{e,pp,qq} = zeros(3,4);
            A_dt{e,pp,qq}  = zeros(4,4);
            A_w{e}         = zeros(4,12);          
            
            % Ciclo para nodos principales
            for i = 1:4
                % Se ensambla la matriz de funciones de forma N
                N{e,pp,qq}(:,[3*i-2 3*i-1 3*i]) = [ NNforma(i)    0             0
                                                    0            NNforma(i)    0
                                                    0            0           NNforma(i) ];
                
                % Se define la matriz B de flexion en los nodos
                dNi_dx = (+dy_deta_N*ddN_dxi(i) - dy_dxi_N*ddN_deta(i))/det_JeN(pp,qq);
                dNi_dy = (-dx_deta_N*ddN_dxi(i) + dx_dxi_N*ddN_deta(i))/det_JeN(pp,qq);
                Bb_t{e,pp,qq}(:,[3*i-2 3*i-1 3*i]) = [  0    dNi_dx    0
                                                        0    0         dNi_dy
                                                        0    dNi_dy    dNi_dx    ];
            end
            
            % Ciclo para nodos secundarios
            for k = 1:4
                % Se definen las distancias de cada lado
                xji(k) = xee(k+1)-xee(k);    yji(k) = yee(k+1)-yee(k);
                Lk(k) = sqrt(xji(k)^2+yji(k)^2);                
                
                % Se definen los senos y cosenos directores
                Ck(k) = xji(k)/Lk(k);  Sk(k) = yji(k)/Lk(k);
                
                % Se define la matriz B de flexion en los nodos secundarios
                dPk_dx = (+dy_deta_N*ddP_dxi(k) - dy_dxi_N*ddP_deta(k))/det_JeN(pp,qq);
                dPk_dy = (-dx_deta_N*ddP_dxi(k) + dx_dxi_N*ddP_deta(k))/det_JeN(pp,qq);
                Bb_dt{e,pp,qq}(:,k) = [ dPk_dx*Ck(k)
                                        dPk_dy*Sk(k)
                                        dPk_dy*Ck(k) + dPk_dx*Sk(k)   ];
            end
            
            % Se define la matriz An
            phi_k = (2/((5/6)*(1-nu))).*(t^2./Lk.^2);
            A_dt{e,pp,qq} = diag(2/3*Lk.*(1+phi_k));
            A_w{e} = [  1    -xji(1)/2    -yji(1)/2    -1    -xji(1)/2    -yji(1)/2    0    0            0            0    0            0
                        0    0            0            1    -xji(2)/2    -yji(2)/2    -1    -xji(2)/2    -yji(2)/2    0    0            0
                        0    0            0            0    0            0            1    -xji(3)/2    -yji(3)/2    -1    -xji(3)/2    -yji(3)/2
                       -1    -xji(4)/2    -yji(4)/2    0    0            0            0    0            0            1    -xji(4)/2    -yji(4)/2];
            A_n = A_dt{e,pp,qq}\A_w{e};
            
            % se ensambla la matriz B de flexion
            Bb{e,pp,qq} = Bb_t{e,pp,qq} + Bb_dt{e,pp,qq}*A_n;
            
            % Se define la matriz B de cortante
            L_phi_5 = Lk(1)*phi_k(1);    L_phi_6 = Lk(2)*phi_k(2);
            L_phi_7 = Lk(3)*phi_k(3);    L_phi_8 = Lk(4)*phi_k(4);
            j_11 = JeN(1,1);            j_12 = JeN(1,2);
            j_21 = JeN(2,1);            j_22 = JeN(2,2);
            Bs_dt{e} = 1/6*[-j_11*(1-eta_gl)*L_phi_5    -j_12*(1+xi_gl)*L_phi_6        j_11*(1+eta_gl)*L_phi_7        j_12*(1-xi_gl)*L_phi_8
                            -j_21*(1-eta_gl)*L_phi_5    -j_22*(1+xi_gl)*L_phi_6        j_21*(1+eta_gl)*L_phi_7        j_22*(1-xi_gl)*L_phi_8    ];
            Bs{e,pp,qq} = Bs_dt{e}*A_n;
            
            % se arma la matriz de rigidez del elemento e por flexion
            Kbe = Kbe + Bb{e,pp,qq}'*Dbe*Bb{e,pp,qq}*det_JeN(pp,qq)*w_gl(pp)*w_gl(qq);
            
            % se arma la matriz de rigidez del elemento e por cortante
            Kse = Kse + Bs{e}'*Dse*Bs{e}*det_JeN(pp,qq)*w_gl(pp)*w_gl(qq);
            
            % vector de fuerzas nodales equivalentes
            if (   (xe(1) >= 0.9999 & xe(1) <= 1.501) & (xe(2) >= 0.9999 & xe(2) <= 1.501) ...
                &  (xe(3) >= 0.9999 & xe(3) <= 1.501) & (xe(4) >= 0.9999 & xe(4) <= 1.501))...
                & ((ye(1) >= 0.9999 & ye(1) <= 2.001) & (ye(2) >= 0.9999 & ye(2) <= 2.001) ...
                &  (ye(3) >= 0.9999 & ye(3) <= 2.001) & (ye(4) >= 0.9999 & ye(4) <= 2.001))
                fe = fe + N{e,pp,qq}'*[q 0 0]'*det_JeN(pp,qq)*w_gl(pp)*w_gl(qq);
            else
                fe = fe + zeros(12,1);
            end
        end
    end
    
    if any(any((det_JeN <= 0) & (det_JeP <= 0)))
        error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
    end
    
    Ke = Kbe + Kse;
    K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
    f(idx{e},:)      = f(idx{e},:)   + fe;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(xnod(:,X) == 0);     lado_y0 = find(xnod(:,Y) == 0);
lado_x2 = find(xnod(:,X) == 2);     lado_y4 = find(xnod(:,Y) == 4);

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
aa = zeros(ngdl,1); aa(c) = ac;  aa(d) = ad; % desplazamientos
qq = zeros(ngdl,1); qq(c) = qd;              % fuerzas nodales equivalentes

vect_mov = reshape(aa,3,nno)'; % vector de movimientos

%% Dibujo la malla de elementos finitos y las deformaciones de esta
escala = 5000; % factor de escalamiento de la deformada
xdef = escala*vect_mov; % posicion de la deformada
figure; 
hold on; 
grid on;
for e = 1:nef
   fill3(xnod(LaG(e,[1 2 3 4 1]),X), ...
         xnod(LaG(e,[1 2 3 4 1]),Y), ...
         xdef(LaG(e,[1 2 3 4 1]),ww),...
         xdef(LaG(e,[1 2 3 4 1]),ww)); %deformada
end
daspect([1 1 1]); % similar a axis equal, pero en 3D
axis tight
colormap jet
title(sprintf('Deformada escalada %d veces',escala),'FontSize',20)
view(3)

%% Se calcula para cada elemento el vector de momentos en los puntos
%% de Gauss
sigma_b = cell(nef,n_gl,n_gl);  % momentos
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            % Se calculan los momentos en los puntos de Gauss
            sigma_b{e,pp,qq} = Dbe*Bb{e,pp,qq}*aa(idx{e});
        end
    end
end

%% Se calcula para cada elemento el vector de cortantes en los puntos
%% de Gauss
QxQy = cell(nef,n_gl,n_gl);  % cortantes
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            % Se calculan las fuerzas cortantes en los puntos de Gauss
            QxQy{e,pp,qq} = Dse*Bs{e}*aa(idx{e});
        end
    end
end

%% Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);

A = [ ... 
    3^(1/2)/2 + 1,          -1/2, 1 - 3^(1/2)/2,          -1/2;
             -1/2, 3^(1/2)/2 + 1,          -1/2, 1 - 3^(1/2)/2;
    1 - 3^(1/2)/2,          -1/2, 3^(1/2)/2 + 1,          -1/2;
             -1/2, 1 - 3^(1/2)/2,          -1/2, 3^(1/2)/2 + 1];

for e = 1:nef                             
   Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + A * [ sigma_b{e,1,1}(1)
											 sigma_b{e,1,2}(1)
											 sigma_b{e,2,1}(1)
											 sigma_b{e,2,2}(1) ];

   My(LaG(e,:),:) = My(LaG(e,:),:)   + A * [ sigma_b{e,1,1}(2)
											 sigma_b{e,1,2}(2)
											 sigma_b{e,2,1}(2)
											 sigma_b{e,2,2}(2) ];
                                        
   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ sigma_b{e,1,1}(3)
											 sigma_b{e,1,2}(3)
											 sigma_b{e,2,1}(3)
											 sigma_b{e,2,2}(3) ];

   Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   + A * [ QxQy{e,1,1}(1)
											 QxQy{e,1,2}(1)
											 QxQy{e,2,1}(1)
											 QxQy{e,2,2}(1) ];
                                            
   Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   + A * [ QxQy{e,1,1}(2)
											 QxQy{e,1,2}(2)
											 QxQy{e,2,1}(2)
											 QxQy{e,2,2}(2) ];
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end 
 
%% Alisado (promedio de los momentos y cortantes en los nodos)
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

%% Se calculan los momentos de disenio de Wood y Armer
[Mxast_sup, Myast_sup, Mxast_inf, Myast_inf] = arrayfun(@WoodArmer, Mx, My, Mxy);
Mmax = max(abs([Mxast_sup; Myast_sup; Mxast_inf; Myast_inf]));

% se graficaran los momentos de disenio utilizando la misma escala de
% colores en valor absoluto, de este modo Mxast_sup=+100 y Mxast_inf=-100 
% tendran el mismo color
figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Mxast_sup,  'Momentos M_x^* sup (N-m/m)');
caxis([0 Mmax]);                                % misma escala de colores
colorbar('ylim', [0 max(Mxast_sup)]);           % rango de colores a mostrar
colormap parula
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_sup,  'Momentos M_y^* sup (N-m/m)');
caxis([0 Mmax]);                                % misma escala de colores
colorbar('ylim', [0 max(Myast_sup)]);           % rango de colores a mostrar
colormap parula

figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Mxast_inf,  'Momentos M_x^* inf (N-m/m)');
colormap parula
caxis([-Mmax 0]);                               % misma escala de colores
oldcmap = colormap; colormap(flipud(oldcmap));  % invierto mapa de colores
colorbar('ylim', [min(Mxast_inf) 0]);           % rango de colores a mostrar
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_inf,  'Momentos M_y^* inf (N-m/m)');
colormap parula
caxis([-Mmax 0]);                               % misma escala de colores
oldcmap = colormap; colormap(flipud(oldcmap));  % invierto mapa de colores
colorbar('ylim', [min(Myast_inf) 0]);           % rango de colores a mostrar

%% Finalmente comparamos los desplazamientos calculados con el MEF y la
%% solucion analitica
u = 0.5; v = 1; xi = 1.25; eta = 1.5;
err = zeros(nno,1);
MEF = zeros(nno,1);
analitica = zeros(nno,1);
for i = 1:nno
   MEF(i) = vect_mov(i,ww);
   analitica(i) = calc_w(xnod(i,X), xnod(i,Y), E, nu, t, 2, 4, q, u, v, xi, eta);
   err(i) = abs((MEF(i)-analitica(i))/analitica(i));
end
disp('Observe que al comparar ambos metodos los errores relativos maximos son')
max(err, [], 'omitnan') % = 0.0027815 =  0.27%
disp('es decir son extremadamente pequenios!!!')

%%
return; % bye, bye!

%%
function plot_M_or_Q(nef, xnod, LaG, variable, texto, angulos)
    X = 1; Y = 2;
    hold on; 
    colorbar;
    for e = 1:nef  
       fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), variable(LaG(e,:)));
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
