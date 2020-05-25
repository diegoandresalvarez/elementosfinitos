% Calculo de los desplazamientos en una placa utilizando la teoria de
% Kirchhoff-Love y el elemento finito de Tocher
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

%% Definimos la geometria de la losa
Mesh_1   % malla no tan refinada
%Mesh_2 % malla muy refinada

% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global

%% Se dibuja la malla de elementos finitos
figure; 
hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1 2 3 1]),X), xnod(LaG(e,[1 2 3 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(LaG(e, :), X));
   cgy(e) = mean(xnod(LaG(e, :), Y));
   text(cgx(e), cgy(e), num2str(e), 'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
title('Malla de elementos finitos');

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)

%% matriz constitutiva
De = (E/(1-nu^2)) * [ 1  nu 0
                      nu 1  0
                      0  0  (1-nu)/2 ];
               
Dbe = (t^3/12)*De;       % matriz constitutiva de flexion generalizada   

D = E*t^3/(12*(1-nu^2)); % rigidez a flexion de la placa   

%% matriz L
L = @(x,y) 	[ 0, 0, 0, 2, 0, 0, 6*x,       2*y,   0 
			  0, 0, 0, 0, 0, 2,   0,       2*x, 6*y 
			  0, 0, 0, 0, 2, 0,   0, 4*x + 4*y,   0];
			  
%% vector p
p = @(x,y)  [ 1; x;	y; x^2; x*y; y^2; x^3; (x^2*y+x*y^2); y^3];
			  
%% se define el orden de la cuadratura
orden = 3;
GP = TriGaussPoints(orden); nGP = size(GP,1);
L2 = GP(:,1);   L3 = GP(:,2);   Wi = GP(:,3)/2;

%% Calculo de Ke y fe
inv_A = cell(nef,1);
for e = 1:nef      % ciclo sobre todos los elementos finitos    
    %% Calculo de la matriz de rigidez Ke
	x1 = xnod(LaG(e,1),X);	y1 = xnod(LaG(e,1),Y);
    x2 = xnod(LaG(e,2),X);  y2 = xnod(LaG(e,2),Y);
	x3 = xnod(LaG(e,3),X);	y3 = xnod(LaG(e,3),Y);
	
	Ae =  0.5*det([ 1 x1 y1      %Area del EF e
					1 x2 y2
					1 x3 y3]);               
    if Ae <= 0
        error('Revise las coordenadas locales del EF %d.\n', e);
    end
    
    iint_LT_Db_L_dA = zeros(9);
    for i = 1:nGP
        x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
        y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
        Le = L(x,y);
        iint_LT_Db_L_dA = iint_LT_Db_L_dA + Le'*Dbe*Le*Wi(i);
    end
    % det_J = 2*Ae;
    iint_LT_Db_L_dA = 2*Ae*iint_LT_Db_L_dA;
    
    inv_A{e} = inv( ... 
          [ 1, x1, y1, x1^2, x1*y1, y1^2,   x1^3, x1^2*y1 + x1*y1^2,   y1^3
            0,  1,  0, 2*x1,    y1,    0, 3*x1^2,    y1^2 + 2*x1*y1,      0
            0,  0,  1,    0,    x1, 2*y1,      0,    x1^2 + 2*y1*x1, 3*y1^2
            1, x2, y2, x2^2, x2*y2, y2^2,   x2^3, x2^2*y2 + x2*y2^2,   y2^3
            0,  1,  0, 2*x2,    y2,    0, 3*x2^2,    y2^2 + 2*x2*y2,      0
            0,  0,  1,    0,    x2, 2*y2,      0,    x2^2 + 2*y2*x2, 3*y2^2
            1, x3, y3, x3^2, x3*y3, y3^2,   x3^3, x3^2*y3 + x3*y3^2,   y3^3
            0,  1,  0, 2*x3,    y3,    0, 3*x3^2,    y3^2 + 2*x3*y3,      0
            0,  0,  1,    0,    x3, 2*y3,      0,    x3^2 + 2*y3*x3, 3*y3^2 ]);        

    % Calculo la matriz de rigidez Ke
    Ke = inv_A{e}'*iint_LT_Db_L_dA*inv_A{e};
    
    %% Calculo del vector de fuerzas nodales equivalentes fe
    if        ((x1 >= 0.9999 && x1 <= 1.501) ...
            && (x2 >= 0.9999 && x2 <= 1.501) ...
            && (x3 >= 0.9999 && x3 <= 1.501) ...
            && (y1 >= 0.9999 && y1 <= 2.001) ...
            && (y2 >= 0.9999 && y2 <= 2.001) ...
            && (y3 >= 0.9999 && y3 <= 2.001))
        int_p_dA = 0;
        for i = 1:nGP
            x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
            y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
            int_p_dA = int_p_dA + p(x,y)*Wi(i);
        end
        fe = inv_A{e}'*q*(2*Ae*int_p_dA);
    else
        fe = zeros(9,1);
    end        
  
    % Ensamblo las contribuciones a las matrices globales
    idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:)];
    K(idx,idx) = K(idx,idx) + Ke;
    f(idx,:)   = f(idx,:)   + fe;
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
aa = zeros(ngdl,1);    aa(c) = ac;  aa(d) = ad; % desplazamientos
qq  = zeros(ngdl,1);   qq(c) = qd;              % fuerzas nodales equivalentes

vect_mov = reshape(aa,3,nno)'; % vector de movimientos

%% Dibujo la malla de elementos finitos y las deformaciones de esta
escala = 5000; % factor de escalamiento de la deformada
xdef = escala*vect_mov; % posicion de la deformada

figure; 
hold on; 
grid on;
for e = 1:nef
   fill3(xnod(LaG(e,[1 2 3 1]),X), ...
         xnod(LaG(e,[1 2 3 1]),Y), ...
         xdef(LaG(e,[1 2 3 1]),ww),...
         xdef(LaG(e,[1 2 3 1]),ww));    % deformada
end
daspect([1 1 1]); % similar a axis equal, pero en 3D
axis tight
colormap jet
title(sprintf('Deformada escalada %d veces',escala),'FontSize',20);
view(3);

%% Se calcula para cada elemento el vector de momentos y cortantes en
%% los puntos de Gauss
sigma_b = cell(nef,nGP);   % momentos flectores Mx y My y torsor Mxy
QxQy    = cell(nef,1);     % cortantes Qx y Qy

QQ = [ ...
             16,          20/3, (4*3^(1/2))/3,            -16,           20/3, -(4*3^(1/2))/3,               0, 8/3, 0
 (16*3^(1/2))/3, (4*3^(1/2))/3,             4, (16*3^(1/2))/3, -(4*3^(1/2))/3,              4, -(32*3^(1/2))/3,   0, 8 ];
      
for e = 1:nef
    x1 = xnod(LaG(e,1),X);    y1 = xnod(LaG(e,1),Y);
    x2 = xnod(LaG(e,2),X);    y2 = xnod(LaG(e,2),Y);
    x3 = xnod(LaG(e,3),X);    y3 = xnod(LaG(e,3),Y);

    idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:)];
    
    for i = 1:nGP
        % calculo de los momentos Mx, My y Mxy
        x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
        y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
        sigma_b{e,i} = -Dbe*L(x,y)*inv_A{e}*aa(idx); 
    end
    % calculo de cortantes Qx y Qy: son constantes para todo el EF
    QxQy{e} = -D*QQ*aa(idx);
end

%% Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);

% matriz de extrapolacion
A = [ 9/4, 5/4, -5/4, -5/4
       -9, 5/2,  5/2,    5
       -9, 5/2,    5,  5/2 ];

for e = 1:nef
    Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + A * [sigma_b{e,1}(1)
                                             sigma_b{e,2}(1)
                                             sigma_b{e,3}(1)
                                             sigma_b{e,4}(1) ];

    My(LaG(e,:),:) = My(LaG(e,:),:)   + A * [sigma_b{e,1}(2)
                                             sigma_b{e,2}(2)
                                             sigma_b{e,3}(2)
                                             sigma_b{e,4}(2) ];
                                        
    Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [sigma_b{e,1}(3)
                                             sigma_b{e,2}(3)
                                             sigma_b{e,3}(3)
                                             sigma_b{e,4}(3) ];
                                             
    Qx(LaG(e,:),:) = Qx(LaG(e,:),:) + QxQy{e}(1);
    Qy(LaG(e,:),:) = Qy(LaG(e,:),:) + QxQy{e}(2);
                                        
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

%% Se calculan y grafican los cortantes m??ximos, junto con su angulo de inclinacion
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
