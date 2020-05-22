clear, clc, close all % borro la memoria, la pantalla y las figuras

% Calculo de los desplazamientos en una placa utilizando la teoria de
% Kirchhoff y el elemento finito de placa MZC.

%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

E  = 210e9;       % modulo de elasticidad del solido (Pa) = 200GPa
nu = 0.3;         % coeficiente de Poisson
t  = 0.05;        % espesor de la losa (m)
p  = -10000;      % carga (N/m^2)

% Definimos la geometria de la losa
losa
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
   line(xnod(LaG(e,[1 2 3 4 1]),X), xnod(LaG(e,[1 2 3 4 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X))/2;
   cgy(e) = (xnod(LaG(e,2),Y) + xnod(LaG(e,3),Y))/2;
   h = text(cgx(e), cgy(e), num2str(e)); set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
%axis equal tight
title('Malla de elementos finitos');

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
a_e = zeros(nef,1);  b_e = zeros(nef,1); % a y b  de cada elemento

% matriz constitutiva
De = E/(1-nu^2)* [ 1  nu 0
                   nu 1  0
                   0  0  (1-nu)/2 ];
               
Dfe = (t^3/12)*De; % matriz constitutiva de flexion generalizada   

for e = 1:nef      % ciclo sobre todos los elementos finitos
   % Calculo de la matriz de rigidez del elemento e    
   x1 = xnod(LaG(e,1),X);
   x2 = xnod(LaG(e,2),X);   y2 = xnod(LaG(e,2),Y);
                            y3 = xnod(LaG(e,3),Y);
   
   a_e(e) = (x2-x1)/2;  a = a_e(e);
   b_e(e) = (y3-y2)/2;  b = b_e(e);

   DD = E*t^3/(12*(1-nu^2));   % rigidez a flexion de la placa   
    
   % Calculo la matriz de rigidez (se calculo con el programa
   % c8_func_forma_MZC.m)
   Ke = DD/(a*b)*[ ...
            b^2/a^2 - nu/5 + a^2/b^2 + 7/10,          (2*nu)/5 + b^2/a^2 + 1/10,          (2*nu)/5 + a^2/b^2 + 1/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             b^2/a^2 - nu/10 + 1/10,      a^2/(2*b^2) - (2*nu)/5 - 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         nu/10 + b^2/(2*a^2) - 1/10,         nu/10 + a^2/(2*b^2) - 1/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      b^2/(2*a^2) - (2*nu)/5 - 1/10,             a^2/b^2 - nu/10 + 1/10
                  (2*nu)/5 + b^2/a^2 + 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                 nu,                  nu/10 - b^2/a^2 - 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,              1/10 - b^2/(2*a^2) - nu/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,           b^2/(2*a^2) - (2*nu)/5 - 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0
                  (2*nu)/5 + a^2/b^2 + 1/10,                                 nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,           a^2/(2*b^2) - (2*nu)/5 - 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,              1/10 - a^2/(2*b^2) - nu/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,                  nu/10 - a^2/b^2 - 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15
        nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             nu/10 - b^2/a^2 - 1/10,      a^2/(2*b^2) - (2*nu)/5 - 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,        - (2*nu)/5 - b^2/a^2 - 1/10,          (2*nu)/5 + a^2/b^2 + 1/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      (2*nu)/5 - b^2/(2*a^2) + 1/10,             a^2/b^2 - nu/10 + 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         1/10 - b^2/(2*a^2) - nu/10,         nu/10 + a^2/(2*b^2) - 1/10
                     b^2/a^2 - nu/10 + 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,             - (2*nu)/5 - b^2/a^2 - 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                -nu,           (2*nu)/5 - b^2/(2*a^2) + 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,              nu/10 + b^2/(2*a^2) - 1/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0
              a^2/(2*b^2) - (2*nu)/5 - 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,               (2*nu)/5 + a^2/b^2 + 1/10,                                -nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,                  nu/10 - a^2/b^2 - 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,              1/10 - a^2/(2*b^2) - nu/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15
    7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         1/10 - b^2/(2*a^2) - nu/10,         1/10 - a^2/(2*b^2) - nu/10,     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      (2*nu)/5 - b^2/(2*a^2) + 1/10,             nu/10 - a^2/b^2 - 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,        - (2*nu)/5 - b^2/a^2 - 1/10,        - (2*nu)/5 - a^2/b^2 - 1/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             nu/10 - b^2/a^2 - 1/10,      (2*nu)/5 - a^2/(2*b^2) + 1/10
                 nu/10 + b^2/(2*a^2) - 1/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,           (2*nu)/5 - b^2/(2*a^2) + 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,             - (2*nu)/5 - b^2/a^2 - 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                 nu,                  b^2/a^2 - nu/10 + 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0
                 nu/10 + a^2/(2*b^2) - 1/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,                  a^2/b^2 - nu/10 + 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,             - (2*nu)/5 - a^2/b^2 - 1/10,                                 nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15,           (2*nu)/5 - a^2/(2*b^2) + 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15
        nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10,      b^2/(2*a^2) - (2*nu)/5 - 1/10,             nu/10 - a^2/b^2 - 1/10, 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5,         nu/10 + b^2/(2*a^2) - 1/10,         1/10 - a^2/(2*b^2) - nu/10,     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10,             b^2/a^2 - nu/10 + 1/10,      (2*nu)/5 - a^2/(2*b^2) + 1/10,         b^2/a^2 - nu/5 + a^2/b^2 + 7/10,          (2*nu)/5 + b^2/a^2 + 1/10,        - (2*nu)/5 - a^2/b^2 - 1/10
              b^2/(2*a^2) - (2*nu)/5 - 1/10, (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15,                                  0,              1/10 - b^2/(2*a^2) - nu/10,         b^2/(3*a^2) - nu/15 + 1/15,                                  0,                  nu/10 - b^2/a^2 - 1/10,     nu/15 + (2*b^2)/(3*a^2) - 1/15,                                  0,               (2*nu)/5 + b^2/a^2 + 1/10, (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15,                                -nu
                     a^2/b^2 - nu/10 + 1/10,                                  0,     nu/15 + (2*a^2)/(3*b^2) - 1/15,              nu/10 + a^2/(2*b^2) - 1/10,                                  0,         a^2/(3*b^2) - nu/15 + 1/15,           (2*nu)/5 - a^2/(2*b^2) + 1/10,                                  0, (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15,             - (2*nu)/5 - a^2/b^2 - 1/10,                                -nu, (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15 ];
     
   % Calculo del vector de fuerzas nodales equivalentes del elemento e
   % Fuerzas superficiales
   if (x1 >= 0.9999 && x2 <= 1.501) && (y2 >= 0.9999 && y3 <= 2.001)
      fe = 4*p*a*b*[1/4;  a/12;  b/12;
                    1/4; -a/12;  b/12; 
                    1/4; -a/12; -b/12;
                    1/4;  a/12; -b/12];
   else
      fe = zeros(12,1);
   end;        
  
   % Ensamblo las contribuciones a las matrices globales
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)];
   K(idx,idx) = K(idx,idx) + Ke;
   f(idx,:)   = f(idx,:)   + fe;
end;

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
qq  = zeros(ngdl,1);   qq(c) = qd;              % fuerzas nodales equivalentes

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
axis tight;	colormap jet
title(sprintf('Deformada escalada %d veces',escala),'FontSize',20);
view(3);

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% xi como en la direccion de eta
n_gl = 2;                 % orden de la cuadratura de Gauss-Legendre

% El comando:
[x_gl, w_gl]  = gausslegendre_quad(n_gl);
% calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre

%% Se calcula para cada elemento el vector de momentos en los puntos
%% de Gauss
esf = cell(nef,n_gl,n_gl);  % momentos
for e = 1:nef
   a = a_e(e); b = b_e(e);

   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:) ];
   
	for pp = 1:n_gl
		for qq = 1:n_gl
			xi = x_gl(pp);		eta = x_gl(qq);
   
			% Se calcula matriz Df*B en los puntos de Gauss
			Df_B = DD/4*[...  % = Df*B
				(3*eta*nu*(xi - 1))/b^2 - (3*xi - 3*eta*xi)/a^2,              ((3*xi - 1)*(eta - 1))/a^2,             (nu*(3*eta - 1)*(xi - 1))/b^2,    (3*xi - 3*eta*xi)/a^2 - (3*eta*nu*(xi + 1))/b^2,              ((3*xi + 1)*(eta - 1))/a^2,           -(nu*(3*eta - 1)*(xi + 1))/b^2,  (3*xi + 3*eta*xi)/a^2 + (3*eta*nu*(xi + 1))/b^2,            -((3*xi + 1)*(eta + 1))/a^2,           -(nu*(3*eta + 1)*(xi + 1))/b^2, - (3*xi + 3*eta*xi)/a^2 - (3*eta*nu*(xi - 1))/b^2,            -((3*xi - 1)*(eta + 1))/a^2,             (nu*(3*eta + 1)*(xi - 1))/b^2
			   (3*nu*xi*(eta - 1))/a^2 - (3*eta - 3*eta*xi)/b^2,           (nu*(3*xi - 1)*(eta - 1))/a^2,                ((3*eta - 1)*(xi - 1))/b^2, - (3*eta + 3*eta*xi)/b^2 - (3*nu*xi*(eta - 1))/a^2,           (nu*(3*xi + 1)*(eta - 1))/a^2,              -((3*eta - 1)*(xi + 1))/b^2, (3*eta + 3*eta*xi)/b^2 + (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi + 1)*(eta + 1))/a^2,              -((3*eta + 1)*(xi + 1))/b^2,  (3*eta - 3*eta*xi)/b^2 - (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi - 1)*(eta + 1))/a^2,                ((3*eta + 1)*(xi - 1))/b^2
					 -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),          ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),       -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b),         ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b) ];
			esf{e,pp,qq} = Df_B*aa(idx);     % Calculo el vector de momentos del elem e
		end;
	end;
end;

%% Se calcula para cada elemento el vector de cortantes en los puntos
%% de Gauss
cort = cell(nef,n_gl,n_gl);  % cortantes
for e = 1:nef
   a = a_e(e); b = b_e(e);

   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:) ];
   
	for pp = 1:n_gl
		for qq = 1:n_gl
			xi = x_gl(pp);		eta = x_gl(qq);
   
			% Se calcula matriz R en los puntos de Gauss
			R =	[ - ((3*eta)/4 - 3/4)/a^3 - (3*eta)/(4*a*b^2),  -(3*eta - 3)/(4*a^3), -(3*eta - 1)/(4*a*b^2), ((3*eta)/4 - 3/4)/a^3 + (3*eta)/(4*a*b^2),  -(3*eta - 3)/(4*a^3), (3*eta - 1)/(4*a*b^2), - ((3*eta)/4 + 3/4)/a^3 - (3*eta)/(4*a*b^2),  (3*eta + 3)/(4*a^3), (3*eta + 1)/(4*a*b^2), ((3*eta)/4 + 3/4)/a^3 + (3*eta)/(4*a*b^2),  (3*eta + 3)/(4*a^3), -(3*eta + 1)/(4*a*b^2)
				    - ((3*xi)/4 - 3/4)/b^3 - (3*xi)/(4*a^2*b), -(3*xi - 1)/(4*a^2*b),    -(3*xi - 3)/(4*b^3),   ((3*xi)/4 + 3/4)/b^3 + (3*xi)/(4*a^2*b), -(3*xi + 1)/(4*a^2*b),    (3*xi + 3)/(4*b^3),   - ((3*xi)/4 + 3/4)/b^3 - (3*xi)/(4*a^2*b), (3*xi + 1)/(4*a^2*b),    (3*xi + 3)/(4*b^3),   ((3*xi)/4 - 3/4)/b^3 + (3*xi)/(4*a^2*b), (3*xi - 1)/(4*a^2*b),    -(3*xi - 3)/(4*b^3)];
			cort{e,pp,qq} = -DD*R*aa(idx);     % Calculo el vector de cortantes del elem e
		end;
	end;
end; 

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
   Mx(LaG(e,:),:) = Mx(LaG(e,:),:)   + A * [esf{e,1,1}(1)
											esf{e,1,2}(1)
											esf{e,2,1}(1)
											esf{e,2,2}(1) ];

   My(LaG(e,:),:) = My(LaG(e,:),:)   + A * [esf{e,1,1}(2)
											esf{e,1,2}(2)
											esf{e,2,1}(2)
											esf{e,2,2}(2) ];
                                        
   Mxy(LaG(e,:),:) = Mxy(LaG(e,:),:) + A * [esf{e,1,1}(3)
											esf{e,1,2}(3)
											esf{e,2,1}(3)
											esf{e,2,2}(3) ];

   Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   + A * [cort{e,1,1}(1)
											cort{e,1,2}(1)
											cort{e,2,1}(1)
											cort{e,2,2}(1) ];
                                            
   Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   + A * [cort{e,1,1}(2)
											cort{e,1,2}(2)
											cort{e,2,1}(2)
											cort{e,2,2}(2) ];
                                          
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
norma = .5; % = Mt_max si quiere proporcional
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

%% Finalmente comparamos los desplazamientos calculados con el MEF y la
%% solucion analitica

u = 0.5; v = 1; xi = 1.25; eta = 1.5;
err = zeros(nno,1);
MEF = zeros(nno,1);
analitica = zeros(nno,1);
for i = 1:nno
   MEF(i) = vect_mov(i,ww);
   analitica(i) = calc_w(xnod(i,X), xnod(i,Y), E, nu, t, 2, 4, p, u, v, xi, eta);
   err(i) = abs((MEF(i)-analitica(i))/analitica(i));
end
disp('Observe que al comparar ambos metodos los errores relativos maximos son')
max(err, [], 'omitnan') % = 0.0027815 =  0.27%
disp('es decir son extremadamente pequeños!!!')

%%
return; % bye, bye!
