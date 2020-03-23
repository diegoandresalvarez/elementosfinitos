%% Programa para el calculo de vigas de Timoshenko

clear, clc, close all               % borrar memoria y pantalla

%% defino las constantes y variables
Y = 1; TH = 2;
%filename = 'viga_Uribe_Escamilla_ej_5_5';
filename = 'viga_con_resortes';
archivo_xlsx = fullfile('..', 'ejemplos', [filename '.xlsx']);

%% se lee la posicion de los nodos
T       = readtable(archivo_xlsx, 'Sheet', 'xnod');
idxNODO = T{:,'nodo'};
xnod    = T{idxNODO,'x'};           % posicion de los nodos
L       = diff(xnod);               % longitud de cada EF

nno  = length(xnod);                % numero de nodos
nef  = nno - 1;                     % numero de elementos finitos (EF)
ngdl = 2*nno;                       % numero de grados de libertad
gdl  = [ (1:2:ngdl)' (2:2:ngdl)' ]; % grados de libertad

%% se leen la matriz de conectividad (LaG), el modulo de elasticidad, las 
%  propiedades del material y las cargas
T     = readtable(archivo_xlsx, 'Sheet', 'LaG_EI_q');
idxEF = T{:,'EF'};
LaG   = T{idxEF,{'NL1','NL2'}};  % definicion de EFs con respecto a nodos
E     = T{idxEF,'E'};            % modulo de elasticidad E del EF
I     = T{idxEF,'I'};            % momento de inercia Iz del EF
G     = T{idxEF,'G'};            % momento de cortante G del EF
Aast  = T{idxEF,'Aast'};         % area reducida Aast del EF
qe    = T{idxEF,{'q1e','q2e'}};  % relacion de las cargas distribuidas
qe(isnan(qe)) = 0;               % reemplace los NaNs con ceros

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
f = zeros(ngdl,1);   % vector de fuerzas nodales equivalentes global
for i = 1:length(idxNODO)
   f(gdl(idxNODO(i), dirfp(i))) = fp(i);
end

%% relacion de los resortes
T       = readtable(archivo_xlsx, 'Sheet', 'resortes');
idxNODO = T{:,'nodo'};
tipores = T{:,'tipo'}; % Y=1 (vertical), TH=2 (rotacional)
kres    = T{:,'k'};    % constante del resorte

%% grados de libertad del desplazamiento conocidos y desconocidos
K = sparse(ngdl,ngdl);  % matriz de rigidez global
n_resortes = length(idxNODO);
for i = 1:n_resortes
   idx = gdl(idxNODO(i), tipores(i));
   K(idx,idx) = kres(i);
end

%% VIGA DE TIMOSHENKO:
% Con los programas "K_exacta_viga_T.m" y "f_exacta_carga_trapezoidal_T.m" 
% se calcularon:
%   Ke     = la matriz de rigidez de flexion del elemento e
%   fe     = el vector de fuerzas nodales equivalentes
%   Bb     = la matriz de deformaciones de flexion
%   N      = matriz de funciones de forma
%   dN_dxi = derivada de la matriz de funciones de forma con respecto a xi

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%% equivalentes global para la viga de Euler-Bernoulli
idx = cell(nef);  % grados de libertad del elemento e
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),Y) gdl(LaG(e,1),TH) gdl(LaG(e,2),Y) gdl(LaG(e,2),TH) ];
   Le = L(e);
   
   % Matriz de rigidez de flexion del elemento e
    beta = (12 * E(e)*I(e))/(Le^2 * G(e)*Aast(e));

    Ke = (E(e)*I(e)/((1 + beta)*Le^3)) * [ ...
                                12,          6*Le,   -12,          6*Le
                              6*Le, (4+beta)*Le^2, -6*Le, (2-beta)*Le^2
                             -12,  -6*Le,    12,  -6*Le
                              6*Le, (2-beta)*Le^2, -6*Le, (4+beta)*Le^2 ];

   % vector de fuerzas nodales equivalentes de una carga trapezoidal 
   fe = [ (  Le*(80*E(e)*I(e)*qe(e,1) + 40*E(e)*I(e)*qe(e,2) + 7*G(e)*Aast(e)*Le^2*qe(e,1) + 3*G(e)*Aast(e)*Le^2*qe(e,2)))/(20* G(e)*Aast(e)*Le^2 + 240*E(e)*I(e) )    % = Y1
          (Le^2*(30*E(e)*I(e)*qe(e,1) + 30*E(e)*I(e)*qe(e,2) + 3*G(e)*Aast(e)*Le^2*qe(e,1) + 2*G(e)*Aast(e)*Le^2*qe(e,2)))/(60*(G(e)*Aast(e)*Le^2 +  12*E(e)*I(e)))    % = M1
          (  Le*(40*E(e)*I(e)*qe(e,1) + 80*E(e)*I(e)*qe(e,2) + 3*G(e)*Aast(e)*Le^2*qe(e,1) + 7*G(e)*Aast(e)*Le^2*qe(e,2)))/(20*(G(e)*Aast(e)*Le^2 +  12*E(e)*I(e)))    % = Y2
         -(Le^2*(30*E(e)*I(e)*qe(e,1) + 30*E(e)*I(e)*qe(e,2) + 2*G(e)*Aast(e)*Le^2*qe(e,1) + 3*G(e)*Aast(e)*Le^2*qe(e,2)))/(60*(G(e)*Aast(e)*Le^2 +  12*E(e)*I(e))) ]; % = M2
   
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

ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% calculo de los momentos flectores
%% (se calcula el momento en el centro de cada elemento finito)
% se reserva la memoria
% recuerde que en cada elemento se calculan los momentos en las raices de 
% los polinomios de Legendre de grado dos
xmom  = zeros(1,nef); % posicion donde se calcula momento flector
mom   = zeros(1,nef); % momento flector
xib   = [ 0 ];        % raices del polinom de Legendre de grado 1 (vect. col)

xcor  = zeros(1,nef); % posicion donde se calcula fuerza cortante
cor   = zeros(1,nef); % fuerza cortante
xis   = [ 0 ];        % raiz del polinomio de Legendre de grado 1 (vect. col)

for e = 1:nef
   Le = L(e);    
    
   % lugar donde se calcula el momento flector y la fuerza cortante
   % (centro del EF)
   xmom(:,e) = Le*xib'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   xcor(:,e) = Le*xis'/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
     
   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = a(idx{e});

   % curvatura kappa y momento flector
   Bb = [0 -1 0 1]/Le;            % matriz de deformacion de flexion
   kappa = Bb*ae;                 % curvatura
   mom(:,e) = E(e)*I(e)*kappa;    % momento flector
   
   % IDEA: para mejorar esta estimacion utilice la matriz Bs del EF de 
   % Timoshenko de 3 nodos
   % gamma_xz y fuerza cortante
   Bs = [ -1/Le  (xis-1)/2  1/Le  -(xis+1)/2 ];
   
   gxz = Bs*ae;                   % gamma_xz  
   cor(e) = -Aast(e)*G(e)*gxz;    % fuerza cortante   
end

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
   xx{e} = L(e)*xi/2 + (xnod(LaG(e,1)) + xnod(LaG(e,2)))/2;
   
   % se calcula el desplazamiento al interior del elemento finito
   ww{e} = Nw*ae;
   
   % se calcula el angulo al interior del elemento finito
   tt{e} = atan(Nt*ae);
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
      fprintf('Nodo %3d W = %12.4g kN, Mx = %12.4g kN-m\n', i, q(i,Y), q(i,TH));
   end
end

%% Grafico la solucion analitica y la solucion por el MEF
%% 1) grafico los desplazamientos de la barra
figure(1)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h1t = plot(xx{e}, ww{e}, 'b-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el desplazamiento')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend('Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 2) grafico los angulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2t = plot(xx{e}, tt{e}, 'b-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el giro')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Giro (rad)')               % titulo del eje Y
legend('Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 3) grafico los momentos
figure(2)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h3t = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [mom(e) mom(e)], 'b-'); % grafico solucion por MEF
end
plot(xmom(:), mom(:), 'bx');% grafico solucion por MEF
title({'Solucion con el MEF para el momento flector',...
   '(el momento positivo es aquel que produce traccion en la fibra inferior)'});
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Momento flector (kN-m)')   % titulo del eje Y
legend('Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 4) grafico la fuerza cortante
subplot(2,1,2);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h4t = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [cor(e) cor(e)], 'b-'); % grafico solucion por MEF
end
%for e = 1:nef % ciclo sobre todos los elementos finitos
%   h4t = plot(xcor(:), cor(:), 'r--');   % grafico solucion por MEF
%end
title('Solucion con el MEF para la fuerza cortante');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Fuerza cortante (kN)')      % titulo del eje Y
legend('Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico


%% Comparacion con la solucion exacta (calculada con MAXIMA y el metodo de
%  las funciones de discontinuidad
if strcmp(filename, 'viga_con_resortes') % OJO solo para b=0.1m y h=0.3m
   dir_txt = fullfile('..', 'ejemplos', 'results_viga_con_resortes_T');
   fid = fopen(fullfile(dir_txt, 'x.txt'));   x = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'Vx.txt'));  V = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'Mx.txt'));  M = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'tx.txt'));  t = str2num(fscanf(fid,'%c')); fclose(fid);
   fid = fopen(fullfile(dir_txt, 'vxx.txt')); v = str2num(fscanf(fid,'%c')); fclose(fid);
  
   figure(1)
   subplot(2,1,1);
   hold on;
   h1tEX = plot(x, v, 'r.');
   legend([h1t, h1tEX], 'EF Timoshenko', 'Timoshenko solucion teorica')
   subplot(2,1,2);
   hold on;
   h2tEX = plot(x, t, 'r.');
   legend([h2t, h2tEX], 'EF Timoshenko', 'Timoshenko solucion teorica')

   figure(2)
   subplot(2,1,1);
   hold on;
   h3tEX = plot(x, M, 'r.');
   legend([h3t, h3tEX], 'EF Timoshenko', 'Timoshenko solucion teorica')
   subplot(2,1,2);
   hold on;
   h4tEX = plot(x, V, 'r.');
   legend([h4t, h4tEX], 'EF Timoshenko', 'Timoshenko solucion teorica')
end

%%
return; % bye, bye!
