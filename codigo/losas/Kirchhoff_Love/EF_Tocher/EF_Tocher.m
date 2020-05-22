clear, clc, close all % borro la memoria, la pantalla y las figuras

% Calculo de los desplazamientos en una placa utilizando la teoria de
% Kirchhoff y el elemento finito de Tocher.

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

E  = 210e9;       % modulo de elasticidad del solido (Pa) = 200GPa
nu = 0.3;         % coeficiente de Poisson
t  = 0.05;        % espesor de la losa (m)
p  = -10000;      % carga (N/m^2)

% Definimos la geometria de la losa
Mesh_2
% Mesh_1 malla no tan refinada
% Mesh_2 malla muy refinada

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
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
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
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)

% matriz constitutiva
De = E/(1-nu^2)* [ 1  nu 0
                   nu 1  0
                   0  0  (1-nu)/2 ];
               
Dfe = (t^3/12)*De; % matriz constitutiva de flexion generalizada   

DD = E*t^3/(12*(1-nu^2));   % rigidez a flexion de la placa   

% matriz L
L = @(x,y) 	[ 0, 0, 0, 2, 0, 0, 6*x,       2*y,   0 
			  0, 0, 0, 0, 0, 2,   0,       2*x, 6*y 
			  0, 0, 0, 0, 2, 0,   0, 4*x + 4*y,   0];
			  
% matriz P
P = @(x,y)  [ 1, x,	y, x^2, x*y, y^2, x^3, (x^2*y+x*y^2), y^3]';
			  
% se define el orden de la cuadratura
orden = 3;
GP = TriGaussPoints(orden); num_pun = size(GP,1);
L2 = GP(:,1);   L3 = GP(:,2);   Wi = GP(:,3)/2;

for e = 1:nef      % ciclo sobre todos los elementos finitos
	% Calculo de la matriz de rigidez del elemento e    
	x1 = xnod(LaG(e,1),X);	y1 = xnod(LaG(e,1),Y);
	x2 = xnod(LaG(e,2),X);  y2 = xnod(LaG(e,2),Y);
	x3 = xnod(LaG(e,3),X);	y3 = xnod(LaG(e,3),Y);
	
	Ae =  0.5*det([ 1 x1 y1      %Area del EF e
					1 x2 y2
					1 x3 y3]);               
	if Ae <= 0
		error('Revise las coordenadas locales del EF %d.\n', e);
    end
	
	LDL = zeros(9);
	
	for i = 1:num_pun
        x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
        y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
		Le = L(x,y);
		LDL = LDL + Le'*Dfe*Le*Wi(i);
	end
	LDL = LDL*2*Ae;
	
	A{e} = 	[ 1, x1, y1, x1^2, x1*y1, y1^2,   x1^3, x1^2*y1 + x1*y1^2,   y1^3
			  0,  1,  0, 2*x1,    y1,    0, 3*x1^2,    y1^2 + 2*x1*y1,      0
			  0,  0,  1,    0,    x1, 2*y1,      0,    x1^2 + 2*y1*x1, 3*y1^2
			  1, x2, y2, x2^2, x2*y2, y2^2,   x2^3, x2^2*y2 + x2*y2^2,   y2^3
			  0,  1,  0, 2*x2,    y2,    0, 3*x2^2,    y2^2 + 2*x2*y2,      0
			  0,  0,  1,    0,    x2, 2*y2,      0,    x2^2 + 2*y2*x2, 3*y2^2
			  1, x3, y3, x3^2, x3*y3, y3^2,   x3^3, x3^2*y3 + x3*y3^2,   y3^3
			  0,  1,  0, 2*x3,    y3,    0, 3*x3^2,    y3^2 + 2*x3*y3,      0
			  0,  0,  1,    0,    x3, 2*y3,      0,    x3^2 + 2*y3*x3, 3*y3^2];

	% Calculo la matriz de rigidez
	Ke = (A{e}'\(A{e}'\LDL)')';     % Ke = inv(A{e}')*LDL*inv(A{e});
	
	% Calculo del vector de fuerzas nodales equivalentes del elemento e
	% Fuerzas superficiales
	if ((x1 >= 0.9999 & x1 <= 1.501) & (x2 >= 0.9999 & x2 <= 1.501) & (x3 >= 0.9999 & x3 <= 1.501))...
        & ((y1 >= 0.9999 & y1 <= 2.001) & (y2 >= 0.9999 & y2 <= 2.001) & (y3 >= 0.9999 & y3 <= 2.001))    
        Pe = 0;
		for i = 1:num_pun
            x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
            y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
			Pe = Pe + P(x,y)*Wi(i);
		end
		fe = A{e}'\(p*2*Ae*Pe);
	else
		fe = zeros(9,1);
	end;        
  
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
aa = zeros(ngdl,1);  aa(c) = ac;  aa(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;              % fuerzas nodales equivalentes

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
         xdef(LaG(e,[1 2 3 1]),ww)); %deformada
end
daspect([1 1 1]); % similar a axis equal, pero en 3D
axis tight; colormap jet
title(sprintf('Deformada escalada %d veces',escala),'FontSize',20);
view(3);

%% Se calcula para cada elemento el vector de momentos y cortantes en
%% los puntos de Gauss
esf = cell(nef,num_pun);  % momentos
cort = cell(nef,num_pun);   % cortantes
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes

% Matriz para el calculo de cortantes
% P = [ 1, x,	y, x^2, x*y, y^2, x^3, (x^2*y+x*y^2), y^3];
% R = [diff(P,x,x,x) + diff(P,x,y,y)
%	   diff(P,y,x,x) + diff(P,y,y,y)];
R = [ 0, 0, 0, 0, 0, 0, 6, 2, 0
      0, 0, 0, 0, 0, 0, 0, 2, 6];
	  
for e = 1:nef
	x1 = xnod(LaG(e,1),X);	y1 = xnod(LaG(e,1),Y);
	x2 = xnod(LaG(e,2),X);  y2 = xnod(LaG(e,2),Y);
	x3 = xnod(LaG(e,3),X);	y3 = xnod(LaG(e,3),Y);

    idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:)];
	
	Df_L_A = zeros(3,9);
 	for i = 1:num_pun
        x = (1-L2(i)-L3(i))*x1 + L2(i)*x2 + L3(i)*x3;
        y = (1-L2(i)-L3(i))*y1 + L2(i)*y2 + L3(i)*y3;
       	Df_L_A = (A{e}'\(L(x,y)'*Dfe'))';		% Df_L_A = Dfe*L(x,y)*inv(A{e});
		esf{e,i} = -Df_L_A*aa(idx);				% calculo de los momentos
        
        % calculo de cortantes
        % los cortantes son constantes para todo el elemento finito
        cort{e,i} = -DD*(A{e}'\R')'*aa(idx);	% cort{e,i} = -DD*R*inv(A{e})*aa(idx);
	end;	
end;

%% Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);

P1 = @(x,y) [x.*0+1 x x.*y y];

for e = 1:nef
	x1 = xnod(LaG(e,1),X);	y1 = xnod(LaG(e,1),Y);
	x2 = xnod(LaG(e,2),X);  y2 = xnod(LaG(e,2),Y);
	x3 = xnod(LaG(e,3),X);	y3 = xnod(LaG(e,3),Y);
	
	x = (1-L2-L3)*x1 + L2*x2 + L3*x3;
    y = (1-L2-L3)*y1 + L2*y2 + L3*y3;
	
    X_elem = [x1 x2 x3]';    Y_elem = [y1 y2 y3]';
    
    A1 = [P1(x,y)]; Pt_A = (A1'\P1(X_elem,Y_elem)')';	% matriz para extrapolar
    Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + Pt_A * [esf{e,1}(1)
												esf{e,2}(1)
												esf{e,3}(1)
												esf{e,4}(1) ];

    My(LaG(e,:),:) = My(LaG(e,:),:)   + Pt_A * [esf{e,1}(2)
												esf{e,2}(2)
												esf{e,3}(2)
												esf{e,4}(2) ];
                                        
	Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + Pt_A * [esf{e,1}(3)
												esf{e,2}(3)
												esf{e,3}(3)
												esf{e,4}(3) ];
											 
	Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   + Pt_A * [cort{e,1}(1)
												cort{e,2}(1)
												cort{e,3}(1)
												cort{e,4}(1) ];

    Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   + Pt_A * [cort{e,1}(2)
												cort{e,2}(2)
												cort{e,3}(2)
												cort{e,4}(2) ];
                                        
	num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los momentos y cortantes en los nodos)
Mx  =  Mx./num_elem_ady;  
My  =  My./num_elem_ady;  
Mxy = Mxy./num_elem_ady;   
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady;

%% Se grafican los momentos
figure;
% Dibujo los momentos Mx
subplot(1,3,1); hold on; colorbar;
for e = 1:nef  
   fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), Mx(LaG(e,:)));
end
axis equal tight;	colormap jet
title('Momentos Mx (N-m/m)','FontSize',20);

% Dibujo los momentos My
subplot(1,3,2); hold on; colorbar;
for e = 1:nef  
   fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), My(LaG(e,:)));
end
axis equal tight;	colormap jet
title('Momentos My (N-m/m)','FontSize',20);

% Dibujo los momentos Mxy
subplot(1,3,3); hold on; colorbar;
for e = 1:nef  
   fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), Mxy(LaG(e,:)));
end
axis equal tight;	colormap jet
title('Momentos Mxy (N-m/m)','FontSize',20);

%% Se grafican los cortantes
figure;
% Dibujo los cortantes Qx
subplot(1,2,1); hold on; colorbar;
for e = 1:nef  
   fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), Qx(LaG(e,:)));
end
axis equal tight;	colormap jet
title('Cortantes Qx (N/m)','FontSize',20);

% Dibujo los cortantes Qy
subplot(1,2,2); hold on; colorbar;
for e = 1:nef  
   fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), Qy(LaG(e,:)));
end
axis equal tight;	colormap jet
title('Cortantes Qy (N/m)','FontSize',20);

%% Se calculan y grafican para cada elemento los momentos principales y
%% sus direcciones

Mf1_xy   = (Mx+My)/2 + sqrt(((Mx-My)/2).^2+Mxy.^2); % momento flector maximo
Mf2_xy   = (Mx+My)/2 - sqrt(((Mx-My)/2).^2+Mxy.^2); % momento flector minimo
Mt_max = (Mf1_xy-Mf2_xy)/2;                         % momento torsion maximo
ang  = 0.5*atan2(2*Mxy, Mx-My); % angulo de inclinacion de Mf1_xy

%% Mf1_xy, Mf2_xy, Mt_max
esc = 0.5; % escala para graficar las flechas

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),Mf1_xy(LaG(e,:)))
end;

% Grafique lineas que indican las direcciones principales de Mf1_xy
norma = 1; % = Mf1_xy si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...   % En el nodo grafique una flecha (linea)
   norma.*cos(ang),norma.*sin(ang),... % indicando la direccion principal de Mf1_xy
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(xnod(:,X),xnod(:,Y),...   % la misma flecha ahora en la otra direccion,
   norma.*cos(ang+pi),norma.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;	colormap jet
title('Mf1_{xy} (N-m/m)','FontSize',26); colorbar

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),Mf2_xy(LaG(e,:)))
end;
% Grafique lineas que indiquen direcciones principales de Mf2_xy
norma = 1; % = Mf2_xy si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...             % flecha indicando la direccion
   norma.*cos(ang+pi/2),norma.*sin(ang+pi/2),... % principal de Mf2_xy
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
   norma.*cos(ang-pi/2),norma.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;	colormap jet
title('Mf2_{xy} (N-m/m)','FontSize',26); colorbar

figure;
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),Mt_max(LaG(e,:)))
end;
% Grafique lineas que indiquen direcciones principales de tau_max,
norma = 1; % = Mt_max si quiere proporcional
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
axis equal tight;	colormap jet
title('Mt_{max} (N-m/m)','FontSize',26); colorbar

%% Se calculan y grafican los cortantes máximos

Q_max = sqrt(Qx.^2+Qy.^2);
ang  = atan2(Qy, Qx); % angulo de inclinacion de Mf1_xy

figure
hold on;
for e = 1:nef
   fill(xnod(LaG(e,:),X),xnod(LaG(e,:),Y),Q_max(LaG(e,:)))
end;

% Grafique lineas que indican las direcciones principales de Q_max
norma = 1; % = Q_max si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...   % En el nodo grafique una flecha (linea)
   norma.*cos(ang),norma.*sin(ang),... % indicando la direccion principal de Q_max
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(xnod(:,X),xnod(:,Y),...   % la misma flecha ahora en la otra direccion,
   norma.*cos(ang+pi),norma.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;	colormap jet
title('Q_{max} (N/m)','FontSize',26); colorbar