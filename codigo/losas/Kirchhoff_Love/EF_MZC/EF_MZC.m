%% Calculo de los desplazamientos en una placa utilizando la teoria de
%% Kirchhoff y el elemento finito de placa MZC.
% Por:
% Diego Andres Alvarez Marin
% Sebastian Jaramillo Moreno

clear, clc, %close all % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
global xnod LaG COLOR_RWB
COLOR_RWB = true; % true = escala rojo/blanco/azul; false = jet no compensado

X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

%% se define la losa a calcular
filename = 'losa_rectangular_libro_solidos_efQ4';
%filename = 'losa_rectangular_libro_solidos_efQ4_sin_rot';
%filename = 'losa_rectangular_libro_solidos_efQ4_levantada';
%filename = 'uniforme_efQ4';
archivo_xlsx = fullfile('..', '..', 'ejemplos', [filename '.xlsx']);

%% se leen las coordenadas de los nodos
T       = readtable(archivo_xlsx, 'Sheet', 'xnod');
idxNODO = T{:,'nodo'};
xnod    = T{idxNODO,{'x','y'}}; % = [x,y]
nno     = size(xnod,1); % numero de nodos

%% se lee la matriz de conectividad (LaG) y la carga distribuida fz
T       = readtable(archivo_xlsx, 'Sheet', 'LaG_fz');
idxEF   = T{:,'EF'};
LaG     = T{idxEF,{'NL1','NL2','NL3','NL4'}};
nef     = size(LaG,1);  % numero de EFs
fz      = T{idxEF, 'fz'};
fz(isnan(fz)) = 0;

%% se definen los apoyos y sus desplazamientos
T       = readtable(archivo_xlsx, 'Sheet', 'restric');
idxNODO = T{:,'nodo'};
dirdesp = T{:,'direccion'};
ac      = T{:,'desplazamiento'}; % desplazamientos conocidos en los apoyos

%% grados de libertad del desplazamiento conocidos y desconocidos
ngdl    = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl     = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

n_apoyos = length(idxNODO);
c = zeros(n_apoyos, 1);      % GDL conocidos, correspondientes a los apoyos
for i = 1:n_apoyos
   c(i) = gdl(idxNODO(i), dirdesp(i));
end
if n_apoyos ~= length(unique(c))
    error('Asignaciones de las restricciones en los apoyos repetidas');
end
d = setdiff((1:ngdl)',c);    % GDL desconocidos y libres

%% relacion de cargas puntuales
T = readtable(archivo_xlsx, 'Sheet', 'carga_punt');
idxNODO = T{:,'nodo'};
dirfp   = T{:,'direccion'};
fp      = T{:,'fuerza_puntual'};

%% se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales 
%  equivalentes global "f"
f = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
for i = 1:length(idxNODO)
   f(gdl(idxNODO(i), dirfp(i))) = fp(i);
end

%% se leen algunas variables
T          = readcell(archivo_xlsx, 'Sheet','varios','Range','B1:B9');
E          = T{1}; % modulo de elasticidad E
nu         = T{2}; % coeficiente de Poisson
rho        = T{3}; % densidad del material
g          = T{4}; % aceleracion de la gravedad
t          = T{5}; % espesor de la losa
U_LONG     = T{6}; % unidades de longitud
U_FUERZA   = T{7}; % unidades de fuerza
U_ESFUERZO = T{8}; % unidades de esfuerzo
ESC_W      = T{9}; % factor de escala para los desplazamientos verticales

peso_propio = rho*g*t;             % peso propio por unidad de area

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

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
a_e = zeros(nef,1);  b_e = zeros(nef,1); % a y b de cada elemento (ancho y alto)
idx = cell(nef,1);     % GDL asociados a cada elemento finito

%% matriz constitutiva
De = (E/(1-nu^2)) * [ 1  nu 0
                      nu 1  0
                      0  0  (1-nu)/2 ];
               
Dbe = (t^3/12)*De; % matriz constitutiva de flexion generalizada   

D = E*t^3/(12*(1-nu^2));   % rigidez a flexion de la placa   

%% Calculo de Ke y fe
for e = 1:nef      % ciclo sobre todos los elementos finitos
   % Calculo de la matriz de rigidez del elemento e    
   x1 = xnod(LaG(e,1),X);
   x2 = xnod(LaG(e,2),X);   y2 = xnod(LaG(e,2),Y);
                            y3 = xnod(LaG(e,3),Y);
   
   a = (x2-x1)/2;  a_e(e) = a;
   b = (y3-y2)/2;  b_e(e) = b;
    
   % Calculo la matriz de rigidez Ke
   % Ke se calculo con el programa func_forma_MZC.m
   Ke = D/(a*b)*[ ...
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
   fe = 4*(fz(e) + peso_propio)*a*b*[ 1/4;  a/12;  b/12
                                      1/4; -a/12;  b/12 
                                      1/4; -a/12; -b/12
                                      1/4;  a/12; -b/12 ];
  
   % Ensamblo las contribuciones a las matrices globales
   idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)];
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e})        + fe;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero')

%% se extraen las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

%% se resuelve el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
aa = zeros(ngdl,1); aa(c) = ac;  aa(d) = ad; % desplazamientos
q  = zeros(ngdl,1);  q(c) = qd;              % fuerzas nodales equivalentes

mat_mov = reshape(aa,3,nno)'; % matriz con movimientos
mat_rea = reshape(q,3,nno)';  % matriz con reacciones
%% Dibujo la malla de elementos finitos y las deformaciones de esta
xdef = ESC_W*mat_mov; % posicion de la deformada
figure; 
if COLOR_RWB
    [min_xdef, max_xdef] = bounds(xdef(:,ww));
    colormap(redwhiteblue(min_xdef, max_xdef));
else
    colormap(jet);
end
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
title(sprintf('Deformada escalada %d veces', ESC_W),'FontSize',20)
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);
view(3)

%% Dibujo las reacciones
mat_rea(mat_rea == 0) = NaN;
figure
subplot(1,3,1); 
stem3(xnod(:,X), xnod(:,Y), mat_rea(:,ww), 'filled');
title(['Reacciones Fz [' U_FUERZA ']']);
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);

subplot(1,3,2); stem3(xnod(:,X), xnod(:,Y), mat_rea(:,tx), 'filled');
title(['Reacciones Mx [' U_FUERZA  ' ' U_LONG ']']);
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);

subplot(1,3,3); stem3(xnod(:,X), xnod(:,Y), mat_rea(:,ty), 'filled');
title(['Reacciones My [' U_FUERZA  ' ' U_LONG ']']);
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);

%% Se calcula para cada elemento el vector de momentos en los puntos
%% de Gauss
n_gl = 2;                          % orden de la cuadratura
x_gl = [ -sqrt(1/3); +sqrt(1/3) ]; % raices del polinomio de Legendre
sigma_b = cell(nef,n_gl,n_gl);     % momentos en cada punto de Gauss
for e = 1:nef
    a = a_e(e); b = b_e(e);
   
    for i = 1:n_gl
        for j = 1:n_gl
            xi = x_gl(i);		eta = x_gl(j);
            
            % Se calcula matriz Db*B en los puntos de Gauss
            % Db_Bb se calculo con el programa func_forma_MZC.m
            Db_Bb = D/4*[...  % = Db*Bb
                (3*eta*nu*(xi - 1))/b^2 - (3*xi - 3*eta*xi)/a^2,              ((3*xi - 1)*(eta - 1))/a^2,             (nu*(3*eta - 1)*(xi - 1))/b^2,    (3*xi - 3*eta*xi)/a^2 - (3*eta*nu*(xi + 1))/b^2,              ((3*xi + 1)*(eta - 1))/a^2,           -(nu*(3*eta - 1)*(xi + 1))/b^2,  (3*xi + 3*eta*xi)/a^2 + (3*eta*nu*(xi + 1))/b^2,            -((3*xi + 1)*(eta + 1))/a^2,           -(nu*(3*eta + 1)*(xi + 1))/b^2, - (3*xi + 3*eta*xi)/a^2 - (3*eta*nu*(xi - 1))/b^2,            -((3*xi - 1)*(eta + 1))/a^2,             (nu*(3*eta + 1)*(xi - 1))/b^2
                (3*nu*xi*(eta - 1))/a^2 - (3*eta - 3*eta*xi)/b^2,           (nu*(3*xi - 1)*(eta - 1))/a^2,                ((3*eta - 1)*(xi - 1))/b^2, - (3*eta + 3*eta*xi)/b^2 - (3*nu*xi*(eta - 1))/a^2,           (nu*(3*xi + 1)*(eta - 1))/a^2,              -((3*eta - 1)*(xi + 1))/b^2, (3*eta + 3*eta*xi)/b^2 + (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi + 1)*(eta + 1))/a^2,              -((3*eta + 1)*(xi + 1))/b^2,  (3*eta - 3*eta*xi)/b^2 - (3*nu*xi*(eta + 1))/a^2,         -(nu*(3*xi - 1)*(eta + 1))/a^2,                ((3*eta + 1)*(xi - 1))/b^2
                -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),          ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), -((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b),       -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b), ((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b),         ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b), ((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b), -((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b) ];
            sigma_b{e,i,j} = Db_Bb*aa(idx{e});
        end
    end
end

%% Se calcula para cada elemento el vector de cortantes en el centro del EF
QxQy = cell(nef,1);  % cortantes
xi = 0; eta = 0;     % centro del EF (punto de Gauss)
for e = 1:nef
    a = a_e(e); b = b_e(e);
    
    % QQ se calculo con el programa func_forma_MZC.m
    QQ = [ ...
          -((3*eta)/4 - 3/4)/a^3 - (3*eta)/(4*a*b^2),  -(3*eta - 3)/(4*a^3), -(3*eta - 1)/(4*a*b^2), ((3*eta)/4 - 3/4)/a^3 + (3*eta)/(4*a*b^2),  -(3*eta - 3)/(4*a^3), (3*eta - 1)/(4*a*b^2), - ((3*eta)/4 + 3/4)/a^3 - (3*eta)/(4*a*b^2),  (3*eta + 3)/(4*a^3), (3*eta + 1)/(4*a*b^2), ((3*eta)/4 + 3/4)/a^3 + (3*eta)/(4*a*b^2),  (3*eta + 3)/(4*a^3), -(3*eta + 1)/(4*a*b^2)
            -((3*xi)/4 - 3/4)/b^3 - (3*xi)/(4*a^2*b), -(3*xi - 1)/(4*a^2*b),    -(3*xi - 3)/(4*b^3),   ((3*xi)/4 + 3/4)/b^3 + (3*xi)/(4*a^2*b), -(3*xi + 1)/(4*a^2*b),    (3*xi + 3)/(4*b^3),   - ((3*xi)/4 + 3/4)/b^3 - (3*xi)/(4*a^2*b), (3*xi + 1)/(4*a^2*b),    (3*xi + 3)/(4*b^3),   ((3*xi)/4 - 3/4)/b^3 + (3*xi)/(4*a^2*b), (3*xi - 1)/(4*a^2*b),    -(3*xi - 3)/(4*b^3) ];
    QxQy{e} = -D*QQ*aa(idx{e});
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

   Qx(LaG(e,:),:) = Qx(LaG(e,:),:)   +       QxQy{e}(1);   
   Qy(LaG(e,:),:) = Qy(LaG(e,:),:)   +       QxQy{e}(2);
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end 

%% Alisado (promedio de los momentos y cortantes en los nodos)
Mx  =  Mx./num_elem_ady;  
My  =  My./num_elem_ady;  
Mxy = Mxy./num_elem_ady;   
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady;  

%% Se grafican los momentos Mx, My, Mxy
unitsM = [U_FUERZA '-' U_LONG '/' U_LONG];
plot_M_or_Q({ Mx,  ['Momentos Mx (' unitsM ')']
              My,  ['Momentos My (' unitsM ')']
              Mxy, ['Momentos Mxy (' unitsM ')'] })

%% Se calculan y grafican para cada elemento los momentos principales y
%% sus direcciones
Mt_max = sqrt(((Mx-My)/2).^2 + Mxy.^2); % momento torsion maximo
Mf1_xy = (Mx+My)/2 + Mt_max;            % momento flector maximo
Mf2_xy = (Mx+My)/2 - Mt_max;            % momento flector minimo
ang  = 0.5*atan2(2*Mxy, Mx-My);         % angulo de inclinacion de Mf1_xy

plot_M_or_Q({ Mf1_xy, ['Mf1_{xy} (' unitsM ')'], { ang }
              Mf2_xy, ['Mf2_{xy} (' unitsM ')'], { ang+pi/2 }
              Mt_max, ['Mt_{max} (' unitsM ')'], { ang+pi/4, ang-pi/4 } })
                     
%% Se calculan y grafican los cortantes Qx, Qy y los Qmaximos, junto con 
%% su angulo de inclinacion
Q_max = hypot(Qx, Qy);
ang   = atan2(Qy, Qx);
unitsQ = [U_FUERZA '/' U_LONG];
plot_M_or_Q({ Qx,    ['Cortantes Qx (' unitsQ ')'],      { }
              Qy,    ['Cortantes Qy (' unitsQ ')'],      { }
              Q_max, ['Cortantes Q_{max} (' unitsQ ')'], { ang } })

%% Se calculan los momentos de disenio de Wood y Armer
[Mxast_sup, Myast_sup, Mxast_inf, Myast_inf] = arrayfun(@WoodArmer, Mx, My, Mxy);
plot_M_or_Q({ Mxast_sup,  ['Momentos M_x^* sup (' unitsM ')']
              Myast_sup,  ['Momentos M_y^* sup (' unitsM ')']
              Mxast_inf,  ['Momentos M_x^* inf (' unitsM ')']
              Myast_inf,  ['Momentos M_y^* inf (' unitsM ')'] } );

%% Se reportan los resultados en un archivo .xlsx
% pandas de python para grabar las tablas es mucho mejor :-\
tabla_aq = array2table([(1:nno)', mat_mov, mat_rea],             ...
    'VariableNames', {'nodo', ['w_' U_LONG], 'tx_rad', 'ty_rad', ...
                      ['q_fz_' U_FUERZA],                        ...
                      ['q_mx_' U_FUERZA '_' U_LONG],             ...
                      ['q_my_' U_FUERZA '_' U_LONG]});

tabla_M = array2table([(1:nno)', Mx, My, Mxy,                              ...
                           Mxast_sup, Myast_sup, Mxast_inf, Myast_inf],    ...
    'VariableNames', {'nodo',                                              ...
                      ['Mx_' U_FUERZA '_' U_LONG '_' U_LONG], 'My', 'Mxy', ...
                      'Mxast_sup', 'Myast_sup', 'Mxast_inf', 'Myast_inf'});

tabla_Q = array2table([(1:nno)', Qx, Qy, Q_max], ...
    'VariableNames', {'nodo', ['Qx_' U_FUERZA '_' U_LONG], 'Qy', 'Q_max'});

warning off;
filename_results = ['resultados_' filename '.xlsx'];
writetable(tabla_aq, filename_results, 'Sheet', 'afq')
writetable(tabla_M,  filename_results, 'Sheet', 'momentos')
writetable(tabla_Q,  filename_results, 'Sheet', 'cortantes')
warning on;

fprintf('Calculo finalizado. En "%s" se guardaron los resultados.\n', filename_results)

%% Finalmente comparamos los desplazamientos calculados con el MEF y la
%% solucion analitica
if strcmp(filename, 'losa_rectangular_libro_solidos_efQ4') || ...
   strcmp(filename, 'losa_rectangular_libro_solidos_efQ4_sin_rot')
        
    u = 0.5; v = 1; xi = 1.25; eta = 1.5;
    qdist = -10;
    err = zeros(nno,1);
    MEF = zeros(nno,1);
    analitica = zeros(nno,1);
    for i = 1:nno
       MEF(i) = mat_mov(i,ww);
       analitica(i) = calc_w(xnod(i,X), xnod(i,Y), E, nu, t, 2, 4, qdist, u, v, xi, eta);
       err(i) = abs((MEF(i)-analitica(i))/analitica(i));
    end
    disp('Observe que al comparar ambos metodos los errores relativos maximos son')
    max(err, [], 'omitnan') % = 0.0027815 =  0.27%
    disp('es decir son extremadamente pequenios!!!')
end

%%
return; % bye, bye!

%%
function plot_M_or_Q(MQ)
    global xnod LaG COLOR_RWB
    X = 1; Y = 2;
    nef = size(LaG, 1);

    [nplots, vars] = size(MQ);
    
    figure
    if COLOR_RWB
        [min_MQ, max_MQ] = bounds(cell2mat(MQ(:,1)));
        colormap(redwhiteblue(min_MQ, max_MQ));
    else
        colormap(jet)
    end
    for i = 1:nplots       
        subplot(1, nplots, i);
        title(MQ{i,2}, 'FontSize',20);
        hold on; 
        for e = 1:nef
           fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), MQ{i,1}(LaG(e,:)), ...
               'EdgeAlpha', 0.05);            % transparencia borde del EF
        end
        axis equal tight
        
        % se ajusta la barra de colores
        if COLOR_RWB
            [min_plot, max_plot] = bounds(MQ{i,1});
            colorbar('ylim', [min_plot, max_plot]); % rango total de colores              
            caxis([min_MQ max_MQ]);                 % rango colores del dibujo
        else
            colorbar
        end

        esc = 0.5;
        if vars == 3
            angulos = MQ{i,3};
            norma = 1; % = MQ{i,1} % si se quiere proporcional
            for j = 1:length(angulos)
                % se indica la flecha de la direccion principal
                quiver(xnod(:,X),xnod(:,Y),...             
                    norma.*cos(angulos{j}), norma.*sin(angulos{j}),... 
                    esc, ...                  % con una escala esc
                    'k',...                   % de color negro
                    'ShowArrowHead','off',... % una flecha sin cabeza
                    'LineWidth',2, ...        % con un ancho de linea 2
                    'Marker','.');            % y en el punto (x,y) poner un punto '.'

                % la misma flecha girada 180 grados
                quiver(xnod(:,X),xnod(:,Y),...             
                    norma.*cos(angulos{j}+pi), norma.*sin(angulos{j}+pi),... 
                    esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.');                    
            end            
        end
    end
end