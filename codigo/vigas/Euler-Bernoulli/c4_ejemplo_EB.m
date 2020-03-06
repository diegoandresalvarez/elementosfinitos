%% Programa para el calculo de vigas de Euler-Bernoulli.
% Viga analizada en el Ejemplo 5-5 de Uribe Escamilla
% Este programa esta particularizado para el EF de longitud 0.1 m

clear, clc, close all         % borrar memoria y pantalla

%% defino las constantes y variables
Y = 1; TH = 2;

L   = 0.1;      % m           % longitud de cada EF
E   = 200e6;    % kPa         % modulo de elasticidad de la barra
b   = 0.3;      % m           % base de la viga
h   = 1.5;      % m           % altura de la viga
A   = b*h;      % m^2         % area transversal de la viga
I   = b*h^3/12; % m^4         % momento de inercia con respecto al eje y

xnod = 0:L:19;                % posicion de los nodos
nno = length(xnod);           % numero de nodos
nef = nno - 1;                % numero de elementos finitos (EF)
ngdl = 2*nno;                 % numero de grados de libertad
gdl = [ (1:2:ngdl)' (2:2:ngdl)' ]; % grados de libertad
LaG = [ (1:(nno-1))' (2:nno)'   ]; % definicion de EFs con respecto a nodos

%% Relacion de cargas distribuidas
p = zeros(nef,1);
p(101:160) = -12; % kN/m     % carga distribuida en [1.0, 1.6] m

%% Relacion de cargas puntuales
f = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
f(gdl( 51,Y)) = -30; % kN   % carga puntual en x = 5 m
f(gdl(191,Y)) = -15; % kN   % carga puntual en x = 19 m

%% VIGA DE EULER-BERNOULLI:
% Con el programa "func_forma_euler_bernoulli.m" se calcularon:
% * Ke     = la matriz de rigidez de flexion del elemento e
% * fe     = el vector de fuerzas nodales equivalentes
% * Bf     = la matriz de deformaciones de flexion
% * N      = matriz de funciones de forma
% * dN_dxi = derivada de la matriz de funciones de forma con respecto a xi

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales 
%% equivalentes global para la viga de Euler-Bernoulli
K = zeros(ngdl);  % matriz de rigidez global
idx = cell(nef);  % grados de libertad del elemento e
EI_L3 = E*I/L^3;    
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),Y) gdl(LaG(e,1),TH) gdl(LaG(e,2),Y) gdl(LaG(e,2),TH) ];
   
   % Matriz de rigidez de flexion del elemento e
   Ke = EI_L3 * [  12,   6*L,  -12,   6*L
                  6*L, 4*L^2, -6*L, 2*L^2
                  -12,  -6*L,   12,  -6*L
                  6*L, 2*L^2, -6*L, 4*L^2 ];

   % vector de fuerzas nodales equivalentes        
   fe = +p(e)*L * [ 1/2; L/12; 1/2; -L/12 ];
   
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

%% Relaciono apoyos
%  gdl           desplazamiento
apoyos = [ ...
   gdl(  1,Y)    0    % m
   gdl(  1,TH)   0    % rad
   gdl(101,Y)    0    % m
   gdl(161,Y)    0 ]; % m  

%% grados de libertad del desplazamiento conocidos y desconocidos
c  = apoyos(:,1);            % GDL conocidos
d =  setdiff((1:ngdl)',c);   % GDL desconocidos
ac = apoyos(:,2);            % desplazamientos conocidos

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

ad = Kdd\(fc - Kdc*ac);      % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% calculo de los momentos flectores
%% (se calcula el momento en el centro de cada elemento finito)
% se reserva la memoria
% recuerde que en cada elemento se calculan los momentos en las raices de 
% los polinomios de Legendre de grado dos
xmom = zeros(2,nef); % posicion donde se calcula
mom  = zeros(2,nef); % momento flector
cor  = zeros(1,nef); % fuerza cortante
%xi = linspace(-1,1,10)'; 
xi = [ -sqrt(1/3); sqrt(1/3) ]; % raices del polinom de Legendre de grado 2

% matriz de deformaciones de flexion
Bbe = [ (6*xi)/L^2, (3*xi - 1)/L, -(6*xi)/L^2, (3*xi + 1)/L ];

dN3_dxi3 = [ 3/2, (3*L)/4, -3/2, (3*L)/4];
for e = 1:nef
   % lugar donde se calcula el momento (centro del EF)
   xmom(:,e) = L*xi'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
     
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});
   
   mom(:,e) = E*I*Bbe*ae;              % momento flector   
   cor(e) = E*I*dN3_dxi3*(8/(L^3))*ae; % fuerza cortante   
end

%% se calculan los desplazamientos al interior de cada EF
nint = 10;           % numero de puntos donde se interpolara dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales

% Matriz de funciones de forma y su derivada
N = [ ...
      xi.^3/4 - (3*xi)/4 + 1/2,                  ...
      -(L*(- xi.^3/4 + xi.^2/4 + xi/4 - 1/4))/2, ...
      - xi.^3/4 + (3*xi)/4 + 1/2,                ...
      -(L*(- xi.^3/4 - xi.^2/4 + xi/4 + 1/4))/2 ];   

dN_dxi = [ ...
      (3*xi.^2)/4 - 3/4,                          ...
      -(L*(- (3*xi.^2)/4 + xi/2 + 1/4))/2,        ...
      3/4 - (3*xi.^2)/4,                          ...   
      (L*((3*xi.^2)/4 + xi/2 - 1/4))/2 ];

xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
ww    = cell(nef,1); % interpol desplazamientos en el elemento
tt    = cell(nef,1); % interpol angulo en el elemento
for e = 1:nef        % ciclo sobre todas los elementos finitos
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});

   % interpola sobre la geometria (coord naturales a geometricas)
   xx{e} = L*xi/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   
   % se calcula el desplazamiento al interior del elemento finito
   ww{e} = N*ae;
   
   % se calcula el angulo al interior del elemento finito
   tt{e} = (dN_dxi*2/L)*ae;
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
title('Solucion con el MEF para el desplazamiento', 'FontSize', 15)
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 2) grafico los angulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2eb = plot(xx{e}, tt{e}, 'b-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el giro', 'FontSize', 15)
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
   '(el momento positivo es aquel que produce traccion en la fibra inferior)'},...
   'FontSize',15);
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
title('Solucion con el MEF para la fuerza cortante', 'FontSize', 15);
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Fuerza cortante (kN)')     % titulo del eje Y
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%%
return; % bye, bye!
