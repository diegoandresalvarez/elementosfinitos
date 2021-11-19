clear, clc, close all        % borro la memoria, la pantalla y las figuras

%% definicion del problema
% Calcule los desplazamientos y las reacciones en el empotramiento 
% de la viga mostrada
% 
% | b (carga distribuida de magnitud b)
% |->->->->->->->->->->->->->->->->
% |====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
% 1    2    3    4          nno-1  nno
% |<----longitud L de la barra---->|   el area transversal de la barra es A

% -----------------------------------------------------------------
% Se usaron tres elementos isoparametricos lagrangianos cuadraticos
% -----------------------------------------------------------------

%% defino las variables
nef = 3;                      % numero de elementos finitos (EF)
nno = 2*nef + 1;              % numero de nodos
ngdl = nno;                   % el # de grados de libertad es el mismo # de nodos
E   = 200e9;    % Pa          % modulo de elasticidad de la barra
A   = (0.01)^2; % m^2         % area transversal de la barra
L   = 2;        % m           % longitud de la barra
b   = 1000;     % N/m         % fuerza axial aplicada sobre cada EF
P   = 250;      % N           % carga nodal al final de la barra

xnod = linspace(0,L,nno);     % posicion de los nodos

le   = repmat(L/nef, nef, 1); % longitud de cada EF

LaG = [1 2 3                  % definicion de EFs con respecto a nodos
       3 4 5
       5 6 7];

%% Parametros de la cuadratura de Gauss-Legendre
n_int_gl = 2;                 % orden de la cuadratura de Gauss-Legendre

% El comando:
[xi_gl,w_gl] = gausslegendre_quad(n_int_gl);
% calcula las raices (xi_gl) y los pesos (w_gl) de polinomios de Legendre

% >> [x_gl,w_gl] = gausslegendre_quad(1)
% x_gl = 0;
% w_gl = 2;
% >> [x_gl,w_gl] = gausslegendre_quad(2)
% x_gl = [  -0.577350269189626;  0.577350269189626 ];
% w_gl = [   1.000000000000000;  1.000000000000000 ];
% >> [x_gl,w_gl] = gausslegendre_quad(3)
% x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
% w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
% >> [x_gl,w_gl] = gausslegendre_quad(4)
% x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
% w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547; 0.347854845137453 ];

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
f(nno) = P;        % relaciono la carga puntual en el nodo "nno"

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = zeros(ngdl);   % matriz de rigidez global
De = E*A;          % matriz constitutiva del elemento 
for e = 1:nef      % ciclo sobre todos los elementos finitos
   idx = LaG(e,:);

   Je = le(e)/2;   % Jacobiano del elemento ( = dx_dxi)
   
   % Calculo las matrices de rigidez y el vector de fuerzas nodales 
   % equivalentes del elemento
   Ke = zeros(3);
   fe = zeros(3,1);
   for m = 1:n_int_gl
      % matriz de deformacion del elemento
      xi = xi_gl(m);
      Be = (1/Je)*[xi-1/2, -2*xi, xi+1/2];
      Ke = Ke + w_gl(m)*Be'*De*Be*Je; % matriz de rigidez del elemento e

      % vector de fuerzas nodales equivalentes
      N = [xi.*(xi-1)/2, (1+xi).*(1-xi), xi.*(xi+1)/2]; % matr. de func. de forma
      fe = fe + w_gl(m)*N'*b*Je; % vector de fuerzas nodales equivalentes
   end

   K(idx,idx) = K(idx,idx) + Ke;
   f(idx,:)   = f(idx,:)   + fe;
end

%% grados de libertad del desplazamiento conocidos y desconocidos
c = 1;    d = 2:nno;

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |  % recuerde que qc=0 (siempre)
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = 0;                 % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% se realizan unos calculos intermedios que necesitaremos mas adelante
nint = 10;           % numero de puntos donde se interpolará dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales
N = [xi.*(xi-1)/2, (1+xi).*(1-xi), xi.*(xi+1)/2]; % matr. de func. de forma
xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
uu    = cell(nef,1); % interpol desplazamientos en el elemento
axial = cell(nef,1); % fuerzas axiales en el elemento
for e = 1:nef        % ciclo sobre todas los elementos finitos
   Je = le(e)/2;     % Jacobiano del elemento ( = dx_dxi)
   Be = (1/Je)*[xi-1/2, -2*xi, xi+1/2]; % matriz de deformacion del elemento

   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = [a(LaG(e,1));    a(LaG(e,2));    a(LaG(e,3))   ]; % = a(LaG(e,:))';
   
   % vector de posiciones de nodos locales
   xe = [xnod(LaG(e,1)); xnod(LaG(e,2)); xnod(LaG(e,3))]; %=xnod(LaG(e,:))';

   xx{e} = N*xe; % interpola sobre la geometría (coord naturales a geométricas)
   uu{e} = N*ae; % interpola sobre los desplazamientos
   
   axial{e} = De*Be*ae; % fuerzas axiales en elemento finito e
end

%% imprimo los resultados
format short g
disp('Desplazamientos (m) = ');                   a
disp('Fuerzas nodales equivalentes (N) OJO! = '); f
disp('Fuerzas nodales de equilibrio (N) = ');     q

%% Grafico la solucion analitica y la solucion por el MEF
%% 1) grafico los desplazamientos de la barra
uexacto = @(x) (-b*x.^2/2 + (P + b*L)*x)/(E*A); % solucion analitica

figure                             % cree un nuevo lienzo
x = linspace(0,L,30);              % 30 puntos unif/ distrib. entre 0 y L
plot(x, uexacto(x), 'rx');         % grafico solucion analitica
hold on;                           % no borre el lienzo 
for e = 1:nef % ciclo sobre todos los elementos finitos
   plot(xx{e}, uu{e}, 'b-'); % grafico solucion por MEF
end
title('Comparacion de la solucion analitica con el MEF para el desplazamiento');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','SouthEast');

%% 2) grafico la carga axial de la barra
Nexacta = @(x) (P + b*(L-x));      % solucion analitica para la carga axial

figure                             % cree un nuevo lienzo
plot(x, Nexacta(x), 'rx');         % grafico solucion analitica
hold on;                           % no borre el lienzo
for e = 1:nef % ciclo sobre todos los elementos finitos
    plot(xx{e}, axial{e}, 'b-'); % grafico solucion por MEF
end
title('Comparacion de la solucion analitica con el MEF para la carga axial');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Carga axial (N)')          % titulo del eje Y
legend('solucion exacta','solucion por el MEF', 'Location','NorthEast');

%%
return; % bye, bye!
