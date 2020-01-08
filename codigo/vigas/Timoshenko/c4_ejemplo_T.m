%% Programa para el calculo de vigas de Timoshenko
% Viga analizada en el Ejemplo 5-5 de Uribe Escamilla

c4_ejemplo_EB;  % cargo el ejemplo de Euler-Bernoulli para comparar

clearvars -except L E b h h?eb % borra variables excepto las mencionadas
clc                            % borrar pantalla

%% defino las constantes y variables
Y = 1; TH = 2;

A   = b*h;      % m^2         % area transversal de la viga
I   = b*h^3/12; % m^4         % momento de inercia con respecto al eje y
Aast = (5/6)*A;               % area reducida = coeficiente de distorsion*A
G   = E/(2*(1+0.3));          % modulo de rigidez

xnod = 0:L:19;                % posicion de los nodos
nno = length(xnod);           % numero de nodos
nef = nno - 1;                % numero de elementos finitos (EF)
ngdl = 2*nno;                 % numero de grados de libertad
gdl = [ (1:2:2*nno)' (2:2:2*nno)' ]; % grados de libertad
LaG = [ (1:(nno-1))' (2:nno)'     ]; % definicion de EFs con respecto a nodos

%% Relacion de cargas distribuidas
p = zeros(nef,1);
p(101:160) = -1200; % N/m     % carga distribuida en [1.0, 1.6] m

%% Relacion de cargas puntuales
f = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
f(gdl( 51,Y)) = -3000; % N   % carga puntual en x = 5 m
f(gdl(191,Y)) = -1500; % N   % carga puntual en x = 18 m

%% VIGA DE TIMOSHENKO:
% Con el programa c4_func_forma_timoshenko_lineal.m se calcularon:
% * Kf = la matriz de rigidez de flexion del elemento e
% * Kc = la matriz de rigidez de cortante del elemento e
% * fe = el vector de fuerzas nodales equivalentes
% * Bf = la matriz de deformaciones de flexion
% * Bc = la matriz de deformaciones de flexion
% * N  = matriz de funciones de forma

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales 
%% equivalentes global para la viga de Timoshenko
K = zeros(ngdl);  % matriz de rigidez global
idx = cell(nef);  % grados de libertad del elemento e
EI_L = E*I/L; GAast_L = G*Aast/L;
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),Y) gdl(LaG(e,1),TH) gdl(LaG(e,2),Y) gdl(LaG(e,2),TH) ];
   
   % Matriz de rigidez de flexion del elemento e
   Kf = EI_L * [ ...
      0  0  0  0
      0 +1  0 -1
      0  0  0  0
      0 -1  0 +1 ];
   
   % Matriz de rigidez de cortante del elemento e   
   Kc = GAast_L * [ ...  % Con cuadratura de GL de orden 1
      1,   L/2    -1    L/2
      L/2  L^2/4  -L/2  L^2/4
     -1   -L/2     1   -L/2
      L/2  L^2/4  -L/2  L^2/4 ];
   
   % Use esta matriz Kc en vez de la anterior si quiere ilustrar el shear
   % locking (bloqueo de la solucion). En este caso calcule la viga con 
   % h = 0.01
%{
   Kc = (G*Aast/L) * [ ...  % Con cuadratura de GL de orden 2
      1    L/2    -1    L/2
     L/2  L^2/3  -L/2  L^2/6
     -1   -L/2     1   -L/2
     L/2  L^2/6  -L/2  L^2/3 ];
%}
    
   Ke = Kf + Kc;

   % vector de fuerzas nodales equivalentes
   fe = +(p(e)*L/2) * [ 1; 0; 1; 0 ];
   
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
end;

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

ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% calculo de los momentos flectores
%% (se calcula el momento en el centro de cada elemento finito)
% se reserva la memoria
% recuerde que en cada elemento se calculan los momentos en las raices de 
% los polinomios de Legendre de grado dos
xmom = zeros(2,nef); % posicion donde se calcula momento flector
mom  = zeros(2,nef); % momento flector
xif = [ -sqrt(1/3); sqrt(1/3) ]; % raices del polinom de Legendre de grado 2

xcor = zeros(2,nef); % posicion donde se calcula fuerza cortante
cor  = zeros(2,nef); % fuerza cortante
xic  = [0; 0];       % raiz del polinomio de Legendre de grado 1

for e = 1:nef
   % lugar donde se calcula el momento (centro del EF)
   xmom(:,e) = L*xif'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   xcor(:,e) = L*xic'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
     
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});

   % curvatura chi y momento flector
   Bf = [0 -1 0 1]/L;          % matriz de deformacion de flexion
   chi = Bf*ae;                % curvatura
   mom(:,e) = E*I*chi;         % momento flector
   
   % gamma_xz y fuerza cortante
   Bc = [ repmat(-1/L,2,1)  (xic-1)/2  repmat(1/L,2,1)  -(xic+1)/2 ];
   
   gxz = Bc*ae;                % gamma_xz  
   cor(:,e) = -Aast*G*gxz;     % fuerza cortante   
end;

%% se calculan los desplazamientos al interior de cada EF
nint = 10;           % numero de puntos donde se interpolara dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales

% Matriz de funciones de forma de desplazamientos y giros
Nw      = [ (1-xi)/2       zeros(nint,1)  (1+xi)/2       zeros(nint,1) ];
Nt      = [ zeros(nint,1)  (1-xi)/2       zeros(nint,1)  (1+xi)/2      ];

xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
ww    = cell(nef,1); % interpol desplazamientos en el elemento
tt    = cell(nef,1); % interpol angulo en el elemento
for e = 1:nef        % ciclo sobre todas los elementos finitos
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});

   % interpola sobre la geometria (coord naturales a geometricas)
   xx{e} = L*xi/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   
   % se calcula el desplazamiento al interior del elemento finito
   ww{e} = Nw*ae;
   
   % se calcula el angulo al interior del elemento finito
   tt{e} = Nt*ae;
end

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(a,2,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: w = %12.4g m, theta = %12.4g rad \n', ...
      i, vect_mov(i,Y), vect_mov(i,TH));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,2,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0])
      fprintf('Nodo %3d W = %12.4g N, Mx = %12.4g N-m\n', i, q(i,Y), q(i,TH));
   end;
end;

%% Grafico la solucion analitica y la solucion por el MEF
%% 1) grafico los desplazamientos de la barra
figure(1)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h1t = plot(xx{e}, ww{e}, 'r--');       % grafico solucion por MEF
end;
title('Solucion con el MEF para el desplazamiento', 'FontSize', 15)
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend([h1eb h1t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');

%% 2) grafico los angulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2t = plot(xx{e}, tt{e}, 'r--');       % grafico solucion por MEF
end;
title('Solucion con el MEF para el giro', 'FontSize', 15)
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Giro (rad)')               % titulo del eje Y
legend([h2eb h2t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');

%% 3) grafico los momentos
figure(2)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo
grid on;                           % reticula
h3t = plot(xmom(:), mom(:), 'r--');% grafico solucion por MEF
title({'Solucion con el MEF para el momento flector',...
   '(el momento positivo es aquel que produce traccion en la fibra inferior)'},...
   'FontSize',15);
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Momento flector (N-m)')    % titulo del eje Y
legend([h3eb h3t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');

%% 4) grafico la fuerza cortante
subplot(2,1,2);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h4t = plot(xcor(:), cor(:), 'r--');   % grafico solucion por MEF
end;
title('Solucion con el MEF para la fuerza cortante', 'FontSize', 15);
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Fuerza cortante (N)')      % titulo del eje Y
legend([h4eb h4t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');

%%
return; % bye, bye!
