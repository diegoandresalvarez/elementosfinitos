% NOTA:
% Tal vez pueda servir de inspiracion:
% https://github.com/JWock82/PyNite/blob/master/PyNite/Quad3D.py

%% 
% Calculo de los desplazamientos en una placa utilizando la teoria de
% Reissner-Mindlin y el elemento finito de placa MITC4
%
% Algoritmo documentado en:
% Bathe, K.-J. and Dvorkin, E.N. (1985), A four-node plate bending element 
% based on Mindlin/Reissner plate theory and a mixed interpolation. Int. J.
% Numer. Meth. Engng., 21: 367-383. https://doi.org/10.1002/nme.1620210213
%
% Por:
% Diego Andres Alvarez Marin (daalvarez@unal.edu.co)

%% borro la memoria, la pantalla y las figuras
clear, clc, %close all 

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

E  = 210e9;       % [Pa]    modulo de elasticidad = 210GPa
nu = 0.3;         %         coeficiente de Poisson
h  = 0.05;        % [m]     espesor de la losa
q  = -10000;      % [N/m^2] carga

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
   
   % Calculo la posicion del centro de gravedad del EF
   cgx(e) = mean(xnod(LaG(e,:),X));
   cgy(e) = mean(xnod(LaG(e,:),Y));
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

%% Se leen las funciones de forma N y P y sus derivadas dN_dr, dN_ds
funciones_de_forma_Bathe_Dvorkin;

%% matrices constitutivas (ec. Bathe-Dvorkin 4)
Db = (E*h^3/(12*(1-nu^2)));   % plate rigidity
Hb = Db * [ 1  nu 0           % matriz constitutiva de flexion generalizada
            nu 1  0           % (Dbe en la nomenclatura del curso) 
            0  0  (1-nu)/2 ]; % (Cb  en la nomenclatura de Bathe-Dvorkin)

G  = E/(2*(1+nu));     % modulo de cortante
Hs = (5/6)*G*h*eye(2); % matriz constitutiva de cortante generalizada (Dse, Cs)

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K   = sparse(ngdl,ngdl);    % matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);        % vector de fuerzas nodales equivalentes global
N   = cell(nef, n_gl, n_gl);
Bb  = cell(nef, n_gl, n_gl);
Bs  = cell(nef, n_gl, n_gl);
idx = cell(nef, 1);         % grados de libertad de cada elemento finito
for e = 1:nef               % ciclo sobre todos los elementos finitos
    %% Longitudes de los lados, cosenos y senos (Figura 4)
    xe = xnod(LaG(e,[3 4 1 2]),X);       ye = xnod(LaG(e,[3 4 1 2]),Y);
    % NOTA: el [3 4 1 2] corresponde al hecho que la numeracion en el 
    % articulo de Bathe y Dvorkin es diferente a aquella con la que se hizo
    % la malla. La malla sigue la numeracion empleada en el articulo de 
    % Katili y es la que se emplea usualmente en elementos finitos
    x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
    x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
    x43 = xe(4) - xe(3);         y43 = ye(4) - ye(3);    
    x14 = xe(1) - xe(4);         y14 = ye(1) - ye(4);
    xji = [ x21 x32 x43 x14 ];   yji = [ y21 y32 y43 y14 ];   
    
    Lk = hypot(xji, yji);      Ck =xji./Lk;      Sk = yji./Lk;
    
    %% Ecuacion Bathe-Dvorkin 10
    Ax = xe(1) - xe(2) - xe(3) + xe(4);
    Bx = xe(1) - xe(2) + xe(3) - xe(4);
    Cx = xe(1) + xe(2) - xe(3) - xe(4);
    Ay = ye(1) - ye(2) - ye(3) + ye(4);
    By = ye(1) - ye(2) + ye(3) - ye(4);
    Cy = ye(1) + ye(2) - ye(3) - ye(4);    
    
    % angulos entre los ejes 
    alpha = atan2(Ay, Ax); % r y x
    beta  = atan2(Cy, Cx); % s y x
    
    %% Ciclo sobre los puntos de Gauss para calcular Kbe, Kse y fe
    Kbe = zeros(12);
    Kse = zeros(12);
    fe  = zeros(12,1);
    det_Je = zeros(n_gl,n_gl); % almacenara los Jacobianos
    
    for pp = 1:n_gl
        for qq = 1:n_gl           
            %% Se evaluan las funciones de forma y sus derivadas en los 
            % puntos de Gauss
            r_gl  = x_gl(pp);            s_gl = x_gl(qq);

            NN     = Nforma(r_gl, s_gl);
            ddN_dr = dN_dr (r_gl, s_gl);       
            ddN_ds = dN_ds (r_gl, s_gl);       
                                   
            %% Matriz jacobiana, su inversa y determinante
            % Se ensambla la matriz jacobiana
            dx_dr = sum(ddN_dr.*xe);   dy_dr = sum(ddN_dr.*ye);
            dx_ds = sum(ddN_ds.*xe);   dy_ds = sum(ddN_ds.*ye);
            
            Je = [ dx_dr   dy_dr
                   dx_ds   dy_ds ];
            
            % Se calcula su inversa
            inv_Je = inv(Je);
            j11 = inv_Je(1,1);              j12 = inv_Je(1,2);
            j21 = inv_Je(2,1);              j22 = inv_Je(2,2);                       
                
            % y su determinante (el Jacobiano)
            det_Je(pp,qq) = det(Je);
            
            %% Se ensambla la matriz de funciones de forma N
            N{e,pp,qq} = zeros(3,12);            
            for i = 1:4
                N{e,pp,qq}(:,[3*i-2 3*i-1 3*i]) = diag([NN(i), NN(i), NN(i)]);
            end            
            
            %% Se calcula la matriz de deformacion por flexion Bb (ec. Katili 40)
            for i = 1:4                
                dNi_dx = j11*ddN_dr(i) + j12*ddN_ds(i); % = ai
                dNi_dy = j21*ddN_dr(i) + j22*ddN_ds(i); % = bi               
                Bb{e,pp,qq}(:,[3*i-2 3*i-1 3*i]) = [ 0   dNi_dx        0
                                                     0        0    dNi_dy
                                                     0   dNi_dy    dNi_dx ];
            end
            
            %% Se calcula la matriz de deformacion por cortante Bs
            % Ecuaciones Bathe-Dvorkin 9
            % La demo de las siguientes ecuaciones esta en "demos_MITC4.m"
            Bgrz = hypot(Cx + r_gl*Bx, Cy + r_gl*By)/(8*det_Je(pp,qq)) *...
                       [        s_gl/2 + 1/2, ...
                          (y21*(s_gl + 1))/4, ...
                         -(x21*(s_gl + 1))/4, ...
                              - s_gl/2 - 1/2, ...
                          (y21*(s_gl + 1))/4, ...
                         -(x21*(s_gl + 1))/4, ...
                                s_gl/2 - 1/2, ...
                          (y43*(s_gl - 1))/4, ...
                         -(x43*(s_gl - 1))/4, ...
                                1/2 - s_gl/2, ...
                          (y43*(s_gl - 1))/4, ...
                         -(x43*(s_gl - 1))/4  ];
            
            Bgsz = hypot(Ax + s_gl*Bx, Ay + s_gl*By)/(8*det_Je(pp,qq)) *...
                       [        r_gl/2 + 1/2, ...
                         -(y14*(r_gl + 1))/4, ...
                          (x14*(r_gl + 1))/4, ...
                                1/2 - r_gl/2, ...
                         -(y32*(r_gl - 1))/4, ...
                          (x32*(r_gl - 1))/4, ...
                                r_gl/2 - 1/2, ...
                         -(y32*(r_gl - 1))/4, ...
                          (x32*(r_gl - 1))/4, ...
                              - r_gl/2 - 1/2, ...
                         -(y14*(r_gl + 1))/4, ...
                          (x14*(r_gl + 1))/4  ];

            % Ecuaciones Bathe-Dvorkin 8
            Bgxz =  Bgrz*sin(beta) - Bgsz*sin(alpha); % 8a
            Bgyz = -Bgrz*cos(beta) + Bgsz*cos(alpha); % 8b
            Bs{e,pp,qq} = [ Bgxz; Bgyz ];

            %% se arma la matriz de rigidez del elemento e por flexion (eq. Katili, 45)
            Kbe = Kbe + Bb{e,pp,qq}'*Hb*Bb{e,pp,qq}*det_Je(pp,qq)*w_gl(pp)*w_gl(qq);
            
            %% se arma la matriz de rigidez del elemento e por cortante (eq. Katili, 47)
            Kse = Kse + Bs{e,pp,qq}'*Hs*Bs{e,pp,qq}*det_Je(pp,qq)*w_gl(pp)*w_gl(qq);
            
            %% vector de fuerzas nodales equivalentes        
            %if (xe(1) >= 0.9999 && xe(2) <= 1.501) && ...
            %   (ye(2) >= 0.9999 && ye(3) <= 2.001)
            if (xe(3) >= 0.9999 && xe(4) <= 1.501) && ...
               (ye(4) >= 0.9999 && ye(1) <= 2.001)            
                fe = fe + N{e,pp,qq}'*[q 0 0]'*det_Je(pp,qq)*w_gl(pp)*w_gl(qq);
            end
        end
    end
    
    %% se verifica que todos los determinantes sean positivos
    if any(det_Je(:) <= 0)
        error('Existen elementos con det(Je(xi,eta)) <= 0 %d.\n', e);
    end
    
    %% ensamblaje matricial
    idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:) ];    
    K(idx{e},idx{e}) = K(idx{e},idx{e}) + Kbe + Kse;
    f(idx{e},:)      = f(idx{e},:)      + fe;
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
%% de Gauss (ecuacion 49)
MxMyMxy = cell(nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            MxMyMxy{e,pp,qq} = Hb*Bb{e,pp,qq}*aa(idx{e});
        end
    end
end

%% Se calcula para cada elemento el vector de cortantes en los puntos
%% de Gauss (ecuacion 50)
QxQy = cell(nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            QxQy{e,pp,qq} = Hs*Bs{e,pp,qq}*aa(idx{e});
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
   3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2
            -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2
   1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1
            -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2 ];

for e = 1:nef                             
    Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + A * [ MxMyMxy{e,1,1}(1)
                                              MxMyMxy{e,1,2}(1)
                                              MxMyMxy{e,2,1}(1)
                                              MxMyMxy{e,2,2}(1) ];

    My(LaG(e,:),:) = My(LaG(e,:),:)   + A * [ MxMyMxy{e,1,1}(2)
                                              MxMyMxy{e,1,2}(2)
                                              MxMyMxy{e,2,1}(2)
                                              MxMyMxy{e,2,2}(2) ];
                                        
    Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [ MxMyMxy{e,1,1}(3)
                                              MxMyMxy{e,1,2}(3)
                                              MxMyMxy{e,2,1}(3)
                                              MxMyMxy{e,2,2}(3) ];
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end 

for e = 1:nef                             

    Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   + A * [ QxQy{e,1,1}(1)
                                              QxQy{e,1,2}(1)
                                              QxQy{e,2,1}(1)
                                              QxQy{e,2,2}(1) ];

    Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   + A * [ QxQy{e,1,1}(2)
                                              QxQy{e,1,2}(2)
                                              QxQy{e,2,1}(2)
                                              QxQy{e,2,2}(2) ];
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

subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_sup,  'Momentos M_y^* sup (N-m/m)');
caxis([0 Mmax]);                                % misma escala de colores
colorbar('ylim', [0 max(Myast_sup)]);           % rango de colores a mostrar

figure
subplot(1,2,1); plot_M_or_Q(nef, xnod, LaG, Mxast_inf,  'Momentos M_x^* inf (N-m/m)');
caxis([-Mmax 0]);                               % misma escala de colores
oldcmap = colormap; colormap(flipud(oldcmap));  % invierto mapa de colores
colorbar('ylim', [min(Mxast_inf) 0]);           % rango de colores a mostrar
subplot(1,2,2); plot_M_or_Q(nef, xnod, LaG, Myast_inf,  'Momentos M_y^* inf (N-m/m)');
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
   analitica(i) = calc_w(xnod(i,X), xnod(i,Y), E, nu, h, 2, 4, q, u, v, xi, eta);
   err(i) = abs((MEF(i)-analitica(i))/analitica(i));
end
disp('Observe que al comparar ambos metodos los errores relativos maximos son')
max(err, [], 'omitnan') % = 0.002128 =  0.21%
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
