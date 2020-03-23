%% Programa para el calculo de vigas de Timoshenko
% Viga analizada en el Ejemplo 5-5 de Uribe Escamilla

c4_ejemplo_EB;  % cargo el ejemplo de Euler-Bernoulli para comparar

clc             % borrar pantalla

%% VIGA DE TIMOSHENKO:
% Con el programa c4_func_forma_timoshenko_lineal.m se calcularon:
% * Kb = la matriz de rigidez de flexion del elemento e
% * Ks = la matriz de rigidez de cortante del elemento e
% * fe = el vector de fuerzas nodales equivalentes
% * Bb = la matriz de deformaciones por flexion
% * Bs = la matriz de deformaciones por cortante
% * N  = matriz de funciones de forma

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales 
%% equivalentes global para la viga de Timoshenko
K = K_ini;
f = f_ini;
for e = 1:nef     % ciclo sobre todos los elementos finitos
   % algunas constantes para hacer el codigo mas legible
   EI_L = E(e)*I(e)/L(e);
   GAast_L = G(e)*Aast(e)/L(e);
   Le = L(e);
   
   % Matriz de rigidez de flexion del elemento e
   Kb = EI_L * [ ...
      0  0  0  0
      0 +1  0 -1
      0  0  0  0
      0 -1  0 +1 ];
   
   % Matriz de rigidez de cortante del elemento e   
   Ks = GAast_L * [ ...  % Con cuadratura de GL de orden 1
      1     Le/2    -1     Le/2
      Le/2  Le^2/4  -Le/2  Le^2/4
     -1    -Le/2     1    -Le/2
      Le/2  Le^2/4  -Le/2  Le^2/4 ];
   
   % Use esta matriz Ks en vez de la anterior si quiere ilustrar el shear
   % locking (bloqueo de la solucion). En este caso calcule la viga con 
   % h = 0.01
%{
   Ks = (G*Aast/L) * [ ...  % Con cuadratura de GL de orden 2
     1      Le/2    -1     Le/2
     Le/2   Le^2/3  -Le/2  Le^2/6
    -1     -Le/2     1    -Le/2
     Le/2   Le^2/6  -Le/2  Le^2/3 ];
%}
    
   Ke = Kb + Ks;

   % vector de fuerzas nodales equivalentes
   fe = +(qe(e,1)*Le/2) * [ 1; 0; 1; 0 ];   % FALTA
   
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
   h1t = plot(xx{e}, ww{e}, 'r-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el desplazamiento')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend([h1eb h1t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 2) grafico los angulos de giro
subplot(2,1,2);
hold on;                           % no borre el lienzo 
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h2t = plot(xx{e}, tt{e}, 'r-');       % grafico solucion por MEF
end
title('Solucion con el MEF para el giro')
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Giro (rad)')               % titulo del eje Y
legend([h2eb h2t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 3) grafico los momentos
figure(2)                          % cree un nuevo lienzo
subplot(2,1,1);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h3t = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [mom(e) mom(e)], 'r-'); % grafico solucion por MEF
end
%h3t = plot(xmom(:), mom(:), 'r--');% grafico solucion por MEF
title({'Solucion con el MEF para el momento flector',...
   '(el momento positivo es aquel que produce traccion en la fibra inferior)'});
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Momento flector (kN-m)')   % titulo del eje Y
legend([h3eb h3t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');
xlim([xnod(1) xnod(end)])          % rango en el eje X del grafico

%% 4) grafico la fuerza cortante
subplot(2,1,2);
hold on;                           % no borre el lienzo
grid on;                           % reticula
for e = 1:nef % ciclo sobre todos los elementos finitos
   h4t = plot([xnod(LaG(e,1)) xnod(LaG(e,2))], [cor(e) cor(e)], 'r-'); % grafico solucion por MEF
end
%for e = 1:nef % ciclo sobre todos los elementos finitos
%   h4t = plot(xcor(:), cor(:), 'r--');   % grafico solucion por MEF
%end
title('Solucion con el MEF para la fuerza cortante');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Fuerza cortante (kN)')      % titulo del eje Y
legend([h4eb h4t], 'Euler-Bernoulli','Timoshenko lineal','Location','Best');
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

   delete([ h1ebEX h2ebEX h3ebEX h4ebEX ]);
   
   figure(1)
   subplot(2,1,1);
   hold on;
   h1tEX = plot(x, v, 'r.');
   legend([h1eb, h1t, h1tEX], 'EF Euler-Bernoulli', 'EF Timoshenko', 'Timoshenko solucion teorica')
   subplot(2,1,2);
   hold on;
   h2tEX = plot(x, t, 'r.');
   legend([h2eb, h2t, h2tEX], 'EF Euler-Bernoulli', 'EF Timoshenko', 'Timoshenko solucion teorica')

   figure(2)
   subplot(2,1,1);
   hold on;
   h3tEX = plot(x, M, 'r.');
   legend([h3eb, h3t, h3tEX], 'EF Euler-Bernoulli', 'EF Timoshenko', 'Timoshenko solucion teorica')
   subplot(2,1,2);
   hold on;
   h4tEX = plot(x, V, 'r.');
   legend([h4eb, h4t, h4tEX], 'EF Euler-Bernoulli', 'EF Timoshenko', 'Timoshenko solucion teorica')
end

%%
return; % bye, bye!
