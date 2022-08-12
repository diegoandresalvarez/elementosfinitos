%% Programa para el calculo de vigas de Euler-Bernoulli.

clear, clc, close all               % borrar memoria y pantalla

%% defino las constantes y variables
Y = 1; TH = 2;
%filename = 'viga_Uribe_Escamilla_ej_5_5';
filename = 'viga_con_resortes';
archivo_xlsx = fullfile('..', 'ejemplos', [filename '.xlsx']);

%% se lee la posicion de los nodos
T       = readtable(archivo_xlsx, 'Sheet', 'xnod', 'Range', 'A:B');
idxNODO = T{:,'nodo'};
xnod    = T{idxNODO,'x'};           % posicion de los nodos
L       = diff(xnod);               % longitud de cada EF

nno  = length(xnod);                % numero de nodos
nef  = nno - 1;                     % numero de elementos finitos (EF)
ngdl = 2*nno;                       % numero de grados de libertad
gdl  = [ (1:2:ngdl)' (2:2:ngdl)' ]; % grados de libertad

%% se leen la matriz de conectividad (LaG), el modulo de elasticidad, las 
%  propiedades del material y las cargas
T     = readtable(archivo_xlsx, 'Sheet', 'LaG_EI_q', 'Range', 'A:I');
idxEF = T{:,'EF'};
LaG   = T{idxEF,{'NL1','NL2'}};  % definicion de EFs con respecto a nodos
E     = T{idxEF,'E'};            % modulo de elasticidad E del EF
I     = T{idxEF,'I'};            % momento de inercia Iz del EF
G     = T{idxEF,'G'};            % modulo de rigidez (para viga de Timoshenko)
Aast  = T{idxEF,'Aast'};         % area de cortante (para viga de Timoshenko)
fz    = T{idxEF,{'q1e','q2e'}};  % relacion de las cargas distribuidas
fz(isnan(fz)) = 0;               % reemplace los NaNs con ceros

%% relacion de los apoyos
T       = readtable(archivo_xlsx, 'Sheet', 'restric');
idxNODO = T{:,'nodo'};
dirdesp = T{:,'direccion'};
ac      = T{:,'desplazamiento'}; % desplazamientos conocidos

%% grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = length(idxNODO);
c = zeros(n_apoyos, 1);          % GDL conocidos
for i = 1:n_apoyos
   c(i,:) = gdl(idxNODO(i), dirdesp(i));
end
d =  setdiff((1:ngdl)',c);       % GDL desconocidos

%% relacion de cargas puntuales
T = readtable(archivo_xlsx, 'Sheet', 'carga_punt');
idxNODO = T{:,'nodo'};
dirfp   = T{:,'direccion'};
fp      = T{:,'fuerza_puntual'};

%% se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales 
%  equivalentes global "f"
f_ini = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
for i = 1:length(idxNODO)
   f_ini(gdl(idxNODO(i), dirfp(i))) = fp(i);
end

%% relacion de los resortes
T       = readtable(archivo_xlsx, 'Sheet', 'resortes');
idxNODO = T{:,'nodo'};
tipores = T{:,'tipo'}; % Y=1 (vertical), TH=2 (rotacional)
kres    = T{:,'k'};    % constante del resorte

%% grados de libertad del desplazamiento conocidos y desconocidos
K_ini = sparse(ngdl,ngdl);  % matriz de rigidez global
n_resortes = length(idxNODO);
for i = 1:n_resortes
   idx = gdl(idxNODO(i), tipores(i));
   K_ini(idx,idx) = kres(i);
end

%% VIGA DE EULER-BERNOULLI:
% Con el programa "func_forma_euler_bernoulli.m" se calcularon:
%   Ke     = la matriz de rigidez de flexion del elemento e
%   fe     = el vector de fuerzas nodales equivalentes
%   Bb     = la matriz de deformaciones de flexion
%   N      = matriz de funciones de forma
%   dN_dxi = derivada de la matriz de funciones de forma con respecto a xi

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%% equivalentes global para la viga de Euler-Bernoulli
idx = cell(nef);  % grados de libertad del elemento e
K = K_ini;
f = f_ini;
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),Y) gdl(LaG(e,1),TH) gdl(LaG(e,2),Y) gdl(LaG(e,2),TH) ];
   Le = L(e);
   
   % Matriz de rigidez de flexion del elemento e
   Ke = (E(e)*I(e)/Le^3) * [  12,   6*Le,   -12,   6*Le
                            6*Le, 4*Le^2, -6*Le, 2*Le^2
                             -12,  -6*Le,    12,  -6*Le
                            6*Le, 2*Le^2, -6*Le, 4*Le^2 ];

   % vector de fuerzas nodales equivalentes de una carga trapezoidal 
   fe = [ (  Le*(7*fz(e,1) + 3*fz(e,2)))/20    % = Y1
          (Le^2*(3*fz(e,1) + 2*fz(e,2)))/60    % = M1
          (  Le*(3*fz(e,1) + 7*fz(e,2)))/20    % = Y2
         -(Le^2*(2*fz(e,1) + 3*fz(e,2)))/60 ]; % = M2
   
   % se ensambla la matriz de rigidez K y el vector de fuerzas nodales
   % equivalentes f
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
end

%% se resuelve el sistema de ecuaciones
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%|   qd   |   | Kcc Kcd || ac |   | fd | 
%|        | = |         ||    | - |    |
%| qc = 0 |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

ad = Kdd\(fc - Kdc*ac);      % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% calculo de los momentos flectores y las fuerzas cortantes
% M = se calcula en las raices del polinomio de GL de orden 2
% V = se calcula en el centro del EF (raiz del polinomio de GL de orden 1)
% se reserva la memoria
xmom = zeros(2,nef); % posicion donde se calcula
mom  = zeros(2,nef); % momento flector
cor  = zeros(1,nef); % fuerza cortante
%xi = linspace(-1,1,10)'; 
xi = [ -sqrt(1/3); sqrt(1/3) ]; % raices del polinom de Legendre de grado 2

for e = 1:nef
   % longitud del elemento finito e
   Le = L(e);
   
   % matriz de deformaciones de flexion
   Bbe = [ (6*xi)/Le^2, (3*xi - 1)/Le, -(6*xi)/Le^2, (3*xi + 1)/Le ];
   
   % lugar donde se calcula el momento (centro del EF)
   xmom(:,e) = Le*xi'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
     
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});
   
   mom(:,e) = E(e)*I(e)*Bbe*ae;                 % momento flector   
   dN3_dxi3 = [ 3/2, (3*Le)/4, -3/2, (3*Le)/4 ];
   cor(e)   = E(e)*I(e)*dN3_dxi3*(8/(Le^3))*ae; % fuerza cortante   
end

%% se calculan los desplazamientos al interior de cada EF
nint = 10;           % numero de puntos donde se interpolara dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales

xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
ww    = cell(nef,1); % interpol desplazamientos en el elemento
tt    = cell(nef,1); % interpol angulo en el elemento
for e = 1:nef        % ciclo sobre todas los elementos finitos
   % longitud del elemento finito e
   Le = L(e);
   
   % Matriz de funciones de forma y su derivada
   N = [ ...
         xi.^3/4 - (3*xi)/4 + 1/2,                   ...
         -(Le*(- xi.^3/4 + xi.^2/4 + xi/4 - 1/4))/2, ...
         - xi.^3/4 + (3*xi)/4 + 1/2,                 ...
         -(Le*(- xi.^3/4 - xi.^2/4 + xi/4 + 1/4))/2 ];   

   dN_dxi = [ ...
         (3*xi.^2)/4 - 3/4,                          ...
         -(Le*(- (3*xi.^2)/4 + xi/2 + 1/4))/2,       ...
         3/4 - (3*xi.^2)/4,                          ...   
         (Le*((3*xi.^2)/4 + xi/2 - 1/4))/2 ];
   
   
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});

   % interpola sobre la geometria (coord naturales a geometricas)
   xx{e} = Le*xi/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   
   % se calcula el desplazamiento al interior del elemento finito
   ww{e} = N*ae;
   
   % se calcula el angulo al interior del elemento finito
   tt{e} = atan((dN_dxi*2/Le)*ae);
end

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(a,2,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: w = %12.4g m, theta = %12.4g rad \n', ...
      i, vect_mov(i,Y), vect_mov(i,TH));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,2,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0])
      fprintf('Nodo %3d Ry = %12.4g kN, Mz = %12.4g kN-m\n', i, q(i,Y), q(i,TH));
   end
end

%% Grafico la solucion analitica y la solucion por el MEF
%% 1) grafico los desplazamientos de la barra
figure(1)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h1eb = plot(xx{e}, ww{e}, 'b-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el desplazamiento')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 2) grafico los angulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2eb = plot(xx{e}, tt{e}, 'b-');% grafico solucion por MEF
end
title('Solucion con el MEF para el giro')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Giro (rad)')               % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 3) grafico los momentos
figure(2)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo
grid on;                           % reticula
h3eb = plot(xmom(:), mom(:), 'b-');% grafico solucion por MEF
title({'Solucion con el MEF para el momento flector',...
   '(el momento positivo es aquel que produce traccion en la fibra inferior)'});
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Momento flector (kN-m)')   % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 4) grafico la fuerza cortante
subplot(2,1,2);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h4eb = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [cor(e) cor(e)], 'b-'); % grafico solucion por MEF
end
title('Solucion con el MEF para la fuerza cortante');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Fuerza cortante (kN)')     % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% Comparacion con la solucion exacta (calculada con MAXIMA y el metodo de
%  las funciones de discontinuidad
if strcmp(filename, 'viga_con_resortes') % OJO solo para b=0.1m y h=0.3m
   dir_txt = fullfile('..', 'ejemplos', 'results_viga_con_resortes_EB');
   fid = fopen(fullfile(dir_txt, 'x.txt'));   x = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'Vx.txt'));  V = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'Mx.txt'));  M = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'tx.txt'));  t = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'vxx.txt')); v = str2num(fscanf(fid,'%c')); fclose(fid);

   figure(1)
   subplot(2,1,1);
   hold on;
   h1ebEX = plot(x, v, 'm.');
   legend([h1eb, h1ebEX], 'Elementos finitos', 'Solucion teorica')
   subplot(2,1,2);
   hold on;
   h2ebEX = plot(x, t, 'm.');
   legend([h2eb, h2ebEX], 'Elementos finitos', 'Solucion teorica')

   figure(2)
   subplot(2,1,1);
   hold on;
   h3ebEX = plot(x, M, 'm.');
   legend([h3eb, h3ebEX], 'Elementos finitos', 'Solucion teorica')
   subplot(2,1,2);
   hold on;
   h4ebEX = plot(x, V, 'm.');
   legend([h4eb, h4ebEX], 'Elementos finitos', 'Solucion teorica')
end

%%
return; % bye, bye!
