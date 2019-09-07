clear, clc, close all        % borro la memoria, la pantalla y las figuras

%% definicion del problema
% Calcule los desplazamientos y las reacciones en el empotramiento 
% de la barra mostrada
% 
% | b (carga distribuida de magnitud b)
% |->->->->->->->->->->->->->->->->
% |====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
% 1    2    3    4          nno-1  nno
% |<----longitud L de la barra---->|   el area transversal de la barra es A

%% defino las variables
nef = 3;                      % numero de elementos finitos (EF)
nno = nef+1;                  % numero de nodos
ngdl = nno;                   % numero de grados de libertad
E   = 200e9;    % Pa          % modulo de elasticidad de la barra
A   = (0.01)^2; % m^2         % area transversal de la barra
L   = 2;        % m           % longitud de la barra
b   = 1000;     % N/m         % fuerza axial aplicada sobre cada EF
P   = 250;      % N           % carga nodal al final de la barra

xnod = linspace(0, L, nno);   % posicion de los nodos
Le  = diff(xnod);             % longitud de cada EF (= repmat(L/nef, nef, 1))
k   = E.*A./Le;               % rigidez de cada EF

LaG = [(1:(nno-1))' (2:nno)']; % definicion de EFs con respecto a nodos

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
f(nno) = P;        % relaciono la carga puntual en el nodo "nno"

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = zeros(ngdl);   % matriz de rigidez global
for e = 1:nef % ciclo sobre todos los elementos finitos
   idx = LaG(e,:);
   K(idx,idx) = K(idx,idx) + k(e)*[1 -1; -1 1];
   f(idx,:)   = f(idx,:)   + ((b*Le(e))/2)*[1; 1];
end;

%% grados de libertad del desplazamiento conocidos y desconocidos
c = 1;    d = 2:ngdl;

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd | 
%|    | = |         ||    | - |    |     Recuerde que qc=0 (siempre)
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = 0;               % desplazamientos conocidos (en el gdl 1)

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% calculo las cargas axiales en cada elemento finito
faxial = zeros(nef,1);
for e = 1:nef % ciclo sobre todas los elementos finitos
   Be = [-1/Le(e) 1/Le(e)];
   ae = [a(LaG(e,1)); a(LaG(e,2))];
   faxial(e) = (E*A)*Be*; % = D*B(e)*a(e)
end;

%% imprimo los resultados
format short g
disp('Desplazamientos (m) = ');                        a 
disp('Fuerzas nodales equivalentes(N) = ');            f
disp('Fuerzas nodales de equilibrio (N) = ');          q
disp('Cargas axiales en cada elemento finito (N) = '); faxial 

%% Grafico la solucion analitica y la solucion por el MEF
% 1) grafico los desplazamientos de la barra
u = @(x) (-b*x.^2/2 + (P + b*L)*x)/(E*A); % solucion analitica para el despl.

figure                             % cree un nuevo lienzo
subplot(2,1,1);                    % grafique en la parte superior (1) del lienzo
xx = linspace(0,L,100);            % 100 puntos unif/ distrib. entre 0 y L
plot(xx, u(xx), 'r');              % grafico solucion analitica
hold on;                           % no borre el lienzo 
plot(xnod, a, 'b.-');              % grafico solucion por MEF
title('Comparacion de la solucion analitica con el MEF para el desplazamiento');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','SouthEast');

% 2) grafico la carga axial de la barra
faxial_exacta = @(x) (P + b*(L-x)); % solucion analitica para la carga axial

subplot(2,1,2);                    % grafique en la parte inferior (2) del lienzo
plot(xx, faxial_exacta(xx), 'r');  % grafico solucion analitica
hold on;                           % no borre el lienzo
for e = 1:nef % ciclo sobre todas los elementos finitos
   plot([xnod(e) xnod(e+1)], [faxial(e) faxial(e)], 'b.-'); % grafico solucion por MEF
end;
title('Comparacion de la solucion analitica con el MEF para la carga axial');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Carga axial (N)')          % titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','NorthEast');

return; %bye, bye!!!
