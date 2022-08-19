%% Programa para el cálculo de vigas de Euler-Bernoulli

%% borrar memoria y pantalla
clear, clc, %close all         

%% se definen las constantes y variables
Y = 1; TH = 2;  % Y: GDL vertical, TH: GDL rotacional
%filename = 'viga_Uribe_Escamilla_ej_5_5';
%filename = 'viga_con_resortes';
filename = 'viga_Winkler';
archivo_xlsx = fullfile('..', 'ejemplos', [filename '.xlsx']);

%% se lee la posición de los nodos (columnas nodo, x)
T       = readtable(archivo_xlsx, 'Sheet', 'xnod', 'Range', 'A:B');
idxNODO = T{:,'nodo'};
xnod    = T{idxNODO,'x'};           % posición de los nodos
L       = diff(xnod);               % longitud de cada EF

nno  = length(xnod);                % número de nodos
nef  = nno - 1;                     % número de elementos finitos (EF)
ngdl = 2*nno;                       % número de grados de libertad
gdl  = [ (1:2:ngdl)' (2:2:ngdl)' ]; % grados de libertad

%% se leen la matriz de conectividad (LaG), el módulo de elasticidad, las 
%  propiedades del material y las cargas
T     = readtable(archivo_xlsx, 'Sheet', 'LaG_EI_q', 'Range', 'A:J');
idxEF = T{:,'EF'};
LaG   = T{idxEF,{'NL1','NL2'}};  % definición de EFs con respecto a nodos
E     = T{idxEF,'E'};            % módulo de elasticidad E del EF
I     = T{idxEF,'I'};            % momento de inercia Iz del EF
G     = T{idxEF,'G'};            % módulo de rigidez (para viga de Timoshenko)
Aast  = T{idxEF,'Aast'};         % área de cortante (para viga de Timoshenko)
fz    = T{idxEF,{'q1e','q2e'}};  % relación de las cargas distribuidas
fz(isnan(fz)) = 0;               % reemplace los NaNs con ceros

%% se leen las constantes de balastro k (cimentacion elastica de Winkler)
bkWinkler = T{idxEF,'bkWinkler'};

%% relación de los apoyos (restricciones)
T       = readtable(archivo_xlsx, 'Sheet', 'restric', 'Range', 'A:C');
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

%% relación de cargas puntuales
T = readtable(archivo_xlsx, 'Sheet', 'carga_punt', 'Range', 'A:C');
idxNODO = T{:,'nodo'};
dirfp   = T{:,'direccion'};
fp      = T{:,'fuerza_puntual'};

%% se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales 
%  equivalentes global "f"
f_ini = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
for i = 1:length(idxNODO)
   f_ini(gdl(idxNODO(i), dirfp(i))) = fp(i);
end

%% relación de los resortes
T       = readtable(archivo_xlsx, 'Sheet', 'resortes', 'Range', 'A:C');
idxNODO = T{:,'nodo'};
tipores = T{:,'tipo'}; % 0=resorte vertical, 1=resorte espiral
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
%   Ke     = la matriz de rigidez de flexión del elemento e
%   fe     = el vector de fuerzas nodales equivalentes
%   Bb     = la matriz de deformaciones de flexión
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

   EI = E(e)*I(e);
   lambda = Le*(bkWinkler(e)/(4*EI))^(1/4);
   den = (sinh(lambda)^2 - sin(lambda)^2);
   nu1 = (4*EI/Le^3)*lambda^3*(sinh(lambda)*cosh(lambda) +  sin(lambda)* cos(lambda))/den;
   nu3 = (2*EI/Le^2)*lambda^2*(sinh(lambda)^2            +  sin(lambda)^2           )/den;
   nu4 = (4*EI/Le^2)*lambda^2*(sinh(lambda)              *  sin(lambda)             )/den;
   nu5 = (2*EI/Le)  *lambda  *(sinh(lambda)*cosh(lambda) -  sin(lambda)* cos(lambda))/den;
   nu6 = (4*EI/Le^3)*lambda^3*(sinh(lambda)* cos(lambda) +  sin(lambda)*cosh(lambda))/den;
   nu7 = (2*EI/Le)  *lambda  *( sin(lambda)*cosh(lambda) - sinh(lambda)* cos(lambda))/den;
   
   Ke = [  nu1  nu3 -nu6  nu4
           nu3  nu5 -nu4  nu7
          -nu6 -nu4  nu1 -nu3
           nu4  nu7 -nu3  nu5 ];

   %{
   % Matriz de rigidez de flexion del elemento e
   Ke = (E(e)*I(e)/Le^3) * [  12,   6*Le,   -12,   6*Le
                            6*Le, 4*Le^2, -6*Le, 2*Le^2
                             -12,  -6*Le,    12,  -6*Le
                            6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
   %}
   % vector de fuerzas nodales equivalentes de una carga trapezoidal 
   fe = [ (  Le*(7*fz(e,1) + 3*fz(e,2)))/20    % = Y1
          (Le^2*(3*fz(e,1) + 2*fz(e,2)))/60    % = M1
          (  Le*(3*fz(e,1) + 7*fz(e,2)))/20    % = Y2
         -(Le^2*(2*fz(e,1) + 3*fz(e,2)))/60 ]; % = M2
   
   % se ensambla la matriz de rigidez K y el vector de fuerzas nodales
   % equivalentes f
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
   %{
   % se suma la contribución de la matriz de rigidez de la cimentación elástica
   if bkWinkler(e) ~= 0
      He = bkWinkler(e)*...
                [   (13*Le)/35, (11*Le^2)/210,      (9*Le)/70, -(13*Le^2)/420
                 (11*Le^2)/210,      Le^3/105,  (13*Le^2)/420,      -Le^3/140
                     (9*Le)/70, (13*Le^2)/420,     (13*Le)/35, -(11*Le^2)/210
                -(13*Le^2)/420,     -Le^3/140, -(11*Le^2)/210,       Le^3/105 ];
      K(idx{e},idx{e}) = K(idx{e},idx{e}) + He;
   end
   %}
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

%% cálculo de los momentos flectores y las fuerzas cortantes
% M = se calcula en las raíces del polinomio de GL de orden 2
% V = se calcula en el centro del EF (raiz del polinomio de GL de orden 1)
% se reserva la memoria
xmom = zeros(2,nef); % posición donde se calcula
mom  = zeros(2,nef); % momento flector  (se calcula en los dos puntos de GL)
cor  = zeros(1,nef); % fuerza cortante  (se calcula en el centro del EF)

% raices del polinomio de Legendre de grado 2
xi = [ -sqrt(1/3); sqrt(1/3) ]; 

for e = 1:nef
   % longitud del elemento finito e
   Le = L(e);
   
   % matriz de deformaciones de flexión
   Bbe = [ (6*xi)/Le^2, (3*xi - 1)/Le, -(6*xi)/Le^2, (3*xi + 1)/Le ];
   
   % lugar donde se calcula el momento (los dos puntos de Gauss del EF)
   xmom(:,e) = Le*xi'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
     
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});
   
   mom(:,e) = E(e)*I(e)*Bbe*ae;                 % momento flector   
   dN3_dxi3 = [ 3/2, (3*Le)/4, -3/2, (3*Le)/4 ];

   % las fuerzas cortantes son constantes dentro del EF (su lectura es más 
   % precisa si se se hace en el centro del EF)
   cor(e)   = E(e)*I(e)*dN3_dxi3*(8/(Le^3))*ae; % fuerza cortante   
end

%% se calculan los desplazamientos al interior de cada EF
nint = 10;           % número de puntos donde se interpolará dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales

xx    = cell(nef,1); % interpolar posiciones (geometría) en el elemento
ww    = cell(nef,1); % interpolar desplazamientos en el elemento
tt    = cell(nef,1); % interpolar ángulo en el elemento
for e = 1:nef        % ciclo sobre todos los elementos finitos
   % longitud del elemento finito e
   Le = L(e);
   
   % Matriz de funciones de forma y sus derivada
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

   % interpola sobre la geometría (coord naturales a geométricas)
   xx{e} = Le*xi/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   
   % se calcula el desplazamiento al interior del elemento finito
   ww{e} = N*ae;
   
   % se calcula el ángulo al interior del elemento finito
   tt{e} = atan((dN_dxi*2/Le)*ae);
end

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                               ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
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

%% Grafico la solución analítica y la solución por el MEF
%% 1) grafico los desplazamientos de la barra
figure(1)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo 
grid on;                           % retícula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h1eb = plot(xx{e}, ww{e}, 'b-');       % grafico solución por el MEF
end
title('Solución con el MEF para el desplazamiento')
xlabel('Eje X (m)')                % título del eje X
ylabel('Desplazamiento (m)')       % título del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del gráfico

%% 2) grafico los ángulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % retícula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2eb = plot(xx{e}, tt{e}, 'b-');% grafico solución por MEF
end
title('Solución con el MEF para el giro')
xlabel('Eje X (m)')                % título del eje X
ylabel('Giro (rad)')               % título del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 3) grafico los momentos
figure(2)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo
grid on;                           % retícula
h3eb = plot(xmom(:), mom(:), 'b-');% grafico solución por MEF
title({'Solución con el MEF para el momento flector',...
   '(el momento positivo es aquel que produce tracción en la fibra inferior)'});
xlabel('Eje X (m)')                % título del eje X
ylabel('Momento flector (kN-m)')   % título del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del gráfico

%% 4) grafico la fuerza cortante
subplot(2,1,2);
hold on;                           % no borre el lienzo
grid on;                           % retícula
for e = 1:nef % ciclo sobre todos los elementos finitos
   % grafico solución por MEF
   h4eb = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [cor(e) cor(e)], 'b-'); 
end
title('Solución con el MEF para la fuerza cortante');
xlabel('Eje X (m)')                % título del eje X
ylabel('Fuerza cortante (kN)')     % título del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del gráfico

%% Comparación con la solución exacta
%  (calculada con MAXIMA y el método de las funciones de discontinuidad
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
   legend([h1eb, h1ebEX], 'Elementos finitos', 'Solución teórica')
   subplot(2,1,2);
   hold on;
   h2ebEX = plot(x, t, 'm.');
   legend([h2eb, h2ebEX], 'Elementos finitos', 'Solución teórica')

   figure(2)
   subplot(2,1,1);
   hold on;
   h3ebEX = plot(x, M, 'm.');
   legend([h3eb, h3ebEX], 'Elementos finitos', 'Solución teórica')
   subplot(2,1,2);
   hold on;
   h4ebEX = plot(x, V, 'm.');
   legend([h4eb, h4ebEX], 'Elementos finitos', 'Solución teórica')
end

%%
return; % bye, bye!
