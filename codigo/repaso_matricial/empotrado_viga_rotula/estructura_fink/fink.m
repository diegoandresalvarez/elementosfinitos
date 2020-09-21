clear, clc, close all

%% constantes
X = 1; Y = 2; TH = 3;

%% Unidades en kN y m
xnod = [ ...  % coordenadas de cada nodo [x, y]
      0      0
    6.0      0
   12.0      0
   18.0      0
    4.5   2.25
    9.0   4.50
   13.5   2.25  ];

% LaG: local a global: matriz que relaciona nodos locales y globales
% fila = barra
% col1 = nodo global asociado a nodo local 1
% col2 = nodo global asociado a nodo local 2
% (se lee la barra x va del nodo i al nodo j)

% nod1  nod2 material  tipo(1=EE,2=ER,3=RR)
barra = [
     2     1     2     2    %  1   2=empotrado, 1=rotula
     2     3     2     1    %  2
     3     4     2     2    %  3
     5     1     1     2    %  4   5=empotrado, 1=rotula
     5     6     1     2    %  5
     7     6     1     2    %  6   7=empotrado, 6=rotula
     7     4     1     2    %  7
     5     2     1     3    %  8
     2     6     1     3    %  9
     6     3     1     3    % 10
     3     7     1     3 ]; % 11
 
LaG = barra(:,[1 2]);  % local a global
mat = barra(:,3);      % material

%    area      inercias_y       modulo de elast.  densidad
%    A(m^2)     I(m^4)          E(kPa)            (kg/m3)
props = [...
    pi*0.04^2  pi*0.04^4/4        200e6          7800     % inclinado
    0.04*0.04  0.04*0.04^3/12     200e6          7800  ]; % horizontal
% seccion circular de radio 4 cm para los elementos inclinados
% seccion rectangular de lado 4 cm para los elementos horizontales
A = props(mat,1); I = props(mat,2); E = props(mat,3); rho = props(mat,4);

%% Se calcula la carga por peso propio asociada a cada barra [kN/m]
wpp = -A.*(rho/1000)*9.81;

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
nbar = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% gdl: grados de libertad
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
% col3 = gdl en direccion angular antihoraria
gdl = [ ...
   1    2  NaN    % 1  % nodos vs grados de libertad
  11   12   13    % 2
  14   15   16    % 3
  17   18  NaN    % 4
   3    4    5    % 5
   6    7  NaN    % 6
   8    9  10 ];  % 7

ngdl = nanmax(gdl(:));  % numero de grados de libertad


%% Se definen los gdl asociados a cada barra
idx = cell(nbar,1);
for e = 1:nbar
   switch barra(e,4) 
      case 1 %EE
         idx{e} = [ gdl(LaG(e,1),[1 2 3]) gdl(LaG(e,2),[1 2 3]) ];    
      case 2 %ER
         idx{e} = [ gdl(LaG(e,1),[1 2 3]) gdl(LaG(e,2),[1 2])   ];             
      case 3 %RR
         idx{e} = [ gdl(LaG(e,1),[1 2])   gdl(LaG(e,2),[1 2])   ];                      
   end   
end

%% cargas aplicadas (gdl carga)
ang = atan2(4.5,3);
cargas_aplica = [ ...
   gdl(2,Y)  -20    
   gdl(3,Y)  -40 
   gdl(5,Y)  -12 
   
   gdl(6,X)  -10*cos(ang)
   gdl(6,Y)  -10*sin(ang)   
   
   gdl(7,X)  -20*cos(ang)
   gdl(7,Y)  -20*sin(ang)   
   
   gdl(4,X)  -10*cos(ang)
   gdl(4,Y)  -10*sin(ang) ];
dofs_cargados  = cargas_aplica(:,1);

fg = zeros(ngdl, 1);
fg(dofs_cargados) = cargas_aplica(:,2);

%% Se calculan las coordenadas x1,y1, x2,y2 de cada barra
x1 = xnod(LaG(:,1), X);  x2 = xnod(LaG(:,2), X);
y1 = xnod(LaG(:,1), Y);  y2 = xnod(LaG(:,2), Y);

%% Se dibuja la estructura junto con su numeracion
figure(1); 
hold on;
cgx = (x1 + x2)/2; cgy = (y1 + y2)/2;   % centro de gravedad
for e = 1:nbar
   line([x1(e) x2(e)], [y1(e) y2(e)]);
   
   % Calculo la posicion del centro de gravedad del triangulo
   h = text(cgx(e), cgy(e), num2str(e)); set(h,'Color', [1 0 0]);
end

axis equal
grid minor
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Numeracion de la estructura');

%% fuerzas nodales equivalentes para las diferentes barras
% (en este ejemplo las fuerzas nodales equivalentes estas siendo 
% especificadas con respecto al sistema de coordenadas globales)
feloc = cell(nbar,1);

%% ensamblo la matriz de rigidez global
Kg   = zeros(ngdl);   % separo memoria
Ke  = cell(nbar,1);
T   = cell(nbar,1);
for e = 1:nbar  % para cada barra   
   % se calcula la matriz de transformacion de coordenadas para la barra e   
   % y la matriz de rigidez local expresada en el sistema de coordenadas 
   % locales para la barra e  
   switch length(idx{e})
      case 4
         [T{e},Kloc] = Te_Keloc_RR(A(e), E(e),       x1(e), y1(e), x2(e), y2(e));
         feloc{e}    = feloc_RR   (wpp(e),           x1(e), y1(e), x2(e), y2(e));
      case 5
         [T{e},Kloc] = Te_Keloc_ER(A(e), E(e), I(e), x1(e), y1(e), x2(e), y2(e));
         feloc{e}    = feloc_ER(wpp(e),              x1(e), y1(e), x2(e), y2(e));
      case 6
         [T{e},Kloc] = Te_Keloc_EE(A(e), E(e), I(e), x1(e), y1(e), x2(e), y2(e));
         feloc{e}    = feloc_EE(wpp(e),              x1(e), y1(e), x2(e), y2(e));         
   end

   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};
   
   Kg(idx{e},idx{e}) = Kg(idx{e},idx{e}) + Ke{e};  % sumo Ke{e} a K global
   fg(idx{e})        = fg(idx{e}) + T{e}'*feloc{e}; % sumo a f global
end

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
apoyos = [...
   gdl(1,Y)  0
   gdl(4,X)  0   
   gdl(4,Y)  0
];

%% introduciendo los soportes inclinados
Tg = speye(ngdl,ngdl);
ang = -30;
Tg(gdl(1,[X Y]),gdl(1,[X Y])) = [ cosd(ang) sind(ang); -sind(ang) cosd(ang) ];
K = Tg*Kg*Tg';

% convierto a sistema de coordenadas con soportes inclinados
f  = Tg*fg;

%%
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

% extraigo las submatrices y especifico las cantidades conocidas
c = apoyos(:,1);
d = setdiff(1:ngdl, c);

Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = apoyos(:,2); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes

%% retorno las fuerzas y los desplazamientos en el sistema de coordenadas
%% donde los grados de libertad son paralelos a los ejes
qg = Tg'*q;
ag = Tg'*a;



%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
%% globales
qe_coord_loc = cell(nbar,1);
for e = 1:nbar % para cada barra
   fprintf('\n\n Fuerzas internas para barra %d en coord. globales = \n', e);
   qe_coord_glob = Ke{e}*ag(idx{e}) - T{e}'*feloc{e};
   disp(qe_coord_glob)
   
   fprintf('\n\n Fuerzas internas para barra %d en coord. locales = \n', e);
   qe_coord_loc{e} = T{e}*qe_coord_glob;
   disp(qe_coord_loc{e});   
end

%% se estiman los movimientos y las reacciones de cada nodo
vect_mov = zeros(nno,3); % vector de movimientos
qq       = zeros(nno,3); % vector de movimientos
for i = 1:nno
   vect_mov(i,1) = ag(gdl(i,1));   qq(i,1) = qg(gdl(i,1));
   vect_mov(i,2) = ag(gdl(i,2));   qq(i,2) = qg(gdl(i,2));
   if isnan(gdl(i,3))
      vect_mov(i,3) = NaN;
      qq(i,3)       = NaN;         
   else
      vect_mov(i,3) = ag(gdl(i,3));      
      qq(i,3)       = qg(gdl(i,3));            
   end
end

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                                                ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g mm, v = %12.4g mm, theta = %12.4g rad \n', ...
      i, 1000*vect_mov(i,X), 1000*vect_mov(i,Y), vect_mov(i,TH));
end


disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)  ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nno   
   if ~isequal(qq(i,:), [0 0 0])
      fprintf('Nodo %3d qx = %12.4g kN, qy = %12.4g, mom = %12.4g kN m\n', ...
         i, qq(i,X), qq(i,Y), qq(i,TH));
   end
end


%% Dibujar la estructura y su deformada
esc_def    = 50;          % escalamiento de la deformada
esc_faxial = 0.05;        % escalamiento del diagrama de axiales
esc_V      = 2;           % escalamiento del diagrama de cortantes
esc_M      = 2;           % escalamiento del diagrama de momentos

xdef = xnod + esc_def*vect_mov(:,[1 2]);

figure(2); hold on; title('Deformada');              xlabel('x, m'); ylabel('y, m'); axis equal
figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal

for e = 1:nbar
   theta = atan2(y2(e)-y1(e),x2(e)-x1(e));
   qxloc_e = @(x) wpp(e)*sin(theta);
   qyloc_e = @(x) wpp(e)*cos(theta);

   switch length(idx{e})
      case 4
      dibujar_barra_deformada_RR(A(e), E(e), I(e), ...
         x1(e),y1(e), x2(e),y2(e), qxloc_e, qyloc_e, qe_coord_loc{e}, T{e}*ag(idx{e}), ...
         esc_def, esc_faxial, esc_V, esc_M);
      case 5
      dibujar_barra_deformada_ER(A(e), E(e), I(e), ...
         x1(e),y1(e), x2(e),y2(e), qxloc_e, qyloc_e, qe_coord_loc{e}, T{e}*ag(idx{e}), ...
         esc_def, esc_faxial, esc_V, esc_M);
      case 6
      dibujar_barra_deformada_EE(A(e), E(e), I(e), ...
         x1(e),y1(e), x2(e),y2(e), qxloc_e, qyloc_e, qe_coord_loc{e}, T{e}*ag(idx{e}), ...
         esc_def, esc_faxial, esc_V, esc_M);
   end
end

%%
return;
