clear, clc, close all

%% Unidades en toneladas y metros
%% constantes
NL1 = 1; NL2 = 2; MAT = 3;
X   = 1; Y   = 2; TH  = 3;

%% Se define la estructura
xnod = [ ...  % coordenadas de cada nodo [x, y]
   3 4
   7 6
   9 0
   0 0 ];

% LaG: local a global: matriz que relaciona nodos locales y globales
% fila = barra
% col1 = nodo global asociado a nodo local 1
% col2 = nodo global asociado a nodo local 2
% (se lee la barra x va del nodo i al nodo j)

barra = [ % NL1  NL2  material
            1    2    1
            4    1    2
            2    3    2 ];

LaG = barra(:, [NL1 NL2]);  % local a global
mat = barra(:, MAT);        % material

%        area      inercias_y       modulo de elasticidad
%        A(m^2)     I(m^4)          E(ton/m^2)
props = [.30*.35   .30*.35^3/12     190e4
         .30*.30   .30*.30^3/12     190e4 ];

A = props(:,1);   I = props(:,2);   E = props(:,3);

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
nbar = size(LaG,1);  % numero de EFs (numero de filas de LaG)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)

%% gdl: grados de libertad
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
% col3 = gdl en direccion angular antihoraria
gdl  = [ (1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)' ]; % nodos vs gdl

%% cargas aplicadas (gdl carga)
cargas_aplica = [ gdl(1,X)  1.5 ];
dofs_cargados = cargas_aplica(:,1);

f = zeros(ngdl, 1);
f(dofs_cargados) = cargas_aplica(:,2);

%% Se dibuja la estructura junto con su numeracion
figure(1); 
hold on;
for e = 1:nbar
   line(xnod(LaG(e,:),X), xnod(LaG(e,:),Y));
   
   % Calculo la posicion del centro de gravedad de la barra
   cgx = (xnod(LaG(e,NL1),X) + xnod(LaG(e,NL2),X))/2;
   cgy = (xnod(LaG(e,NL1),Y) + xnod(LaG(e,NL2),Y))/2;   
   h = text(cgx, cgy, num2str(e)); set(h, 'Color', [1 0 0]);
end

axis equal
grid minor
plot(xnod(:,X), xnod(:,Y), 'ro');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Numeracion de la estructura');

%% fuerzas distribuidas aplicadas sobre las barras en coordenadas locales
ang1 = atan2(2,4);

qxloc = { ...
   @(x) -2.8*sin(ang1)*cos(ang1);
   @(x) 0;
   @(x) 0;
};
qyloc = { ...
   @(x) -2.8*cos(ang1)^2;   
   @(x) 0;
   @(x) 0;
};

%% fuerzas nodales equivalentes para las diferentes barras
% (en este ejemplo las fuerzas nodales equivalentes estas siendo 
% especificadas con respecto al sistema de coordenadas globales)
fe = cell(nbar,1);
for e = 1:nbar
   x1 = xnod(LaG(e,NL1), X);  x2 = xnod(LaG(e,NL2), X);
   y1 = xnod(LaG(e,NL1), Y);  y2 = xnod(LaG(e,NL2), Y);
   fe{e} = calc_fuerzas_nodales_equivalentes(...
        A(mat(e)), E(mat(e)), I(mat(e)), x1,y1, x2,y2, qxloc{e},qyloc{e});
end
%          fxi    fyi    mi     fxj    fyj    mj
%          ton    ton    ton-m  ton    ton    ton-m
% fe{1} = [0     -5.60  -3.733  0     -5.60   +3.733 ]'; % OJO con los signos
% fe{2} = [0      0      0      0      0      0      ]'; % mirar pag 613
% fe{3} = [0      0      0      0      0      0      ]';

%% ensamblo la matriz de rigidez global
K   = zeros(ngdl);   % separo memoria
Ke  = cell(nbar,1);
T   = cell(nbar,1);
idx = zeros(nbar,6);
for e = 1:nbar  % para cada barra
   % saco los 6 gdls de la barra e
   idx(e,:) = [gdl(LaG(e,NL1),:) gdl(LaG(e,NL2),:)];
   
   x1 = xnod(LaG(e,NL1), X);  x2 = xnod(LaG(e,NL2), X);
   y1 = xnod(LaG(e,NL1), Y);  y2 = xnod(LaG(e,NL2), Y);
   
   L = hypot(x2-x1, y2-y1);
   
   % matriz de transformacion de coordenadas para la barra e
   c = (x2-x1)/L;   s = (y2-y1)/L;  % seno y coseno de la inclinacion
   T{e} = [ c  s  0  0  0  0
           -s  c  0  0  0  0
            0  0  1  0  0  0
            0  0  0  c  s  0
            0  0  0 -s  c  0
            0  0  0  0  0  1];
         
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
   AE = A(mat(e))*E(mat(e));       L2=L^2;
   EI = E(mat(e))*I(mat(e));       L3=L^3;
   Kloc = [ AE/L   0         0        -AE/L    0          0       
            0     12*EI/L3   6*EI/L2   0     -12*EI/L3   6*EI/L2
            0      6*EI/L2   4*EI/L    0      -6*EI/L2   2*EI/L
           -AE/L   0         0         AE/L    0         0
            0    -12*EI/L3  -6*EI/L2   0      12*EI/L3  -6*EI/L2
            0      6*EI/L2   2*EI/L    0      -6*EI/L2   4*EI/L];

   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};
   K(idx(e,:),idx(e,:)) = K(idx(e,:),idx(e,:)) + Ke{e}; % ensambla Ke{e} en K global
   f(idx(e,:))          = f(idx(e,:))          + fe{e}; % ensambla fe{e} en f global
end

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
apoyos = [...
   gdl(3,X)  0
   gdl(3,Y)  0
   gdl(3,TH) 0
   gdl(4,X)  0
   gdl(4,Y)  0
   gdl(4,TH) 0
];

c = apoyos(:,1);
d = setdiff(1:ngdl, c);

%%
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = apoyos(:,2); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes


%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
%% globales
qe_loc  = cell(nbar,1);
qe_glob = cell(nbar,1);
for e = 1:nbar % para cada barra
   fprintf('\n\n Fuerzas internas para barra %d en coord. globales = \n', e);
   qe_glob{e} = Ke{e}*a(idx(e,:)) - fe{e};
   disp(qe_glob{e})
   
   fprintf('\n\n Fuerzas internas para barra %d en coord. locales = \n', e);
   qe_loc{e} = T{e}*qe_glob{e};
   disp(qe_loc{e});   
end

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                                            ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
vect_mov = reshape(a,3,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g mm, v = %12.4g mm, theta = %12.4g rad \n', ...
      i, 1000*vect_mov(i,X), 1000*vect_mov(i,Y), vect_mov(i,TH));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo imprimo los diferentes de cero)')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
qq = reshape(q,3,nno)';
for i = 1:nno   
   if ~all(abs(qq(i,:) - [0 0 0]) < 1e-5)
      fprintf('Nodo %3d qx = %12.4g ton, qy = %12.4g ton, mom = %12.4g ton*m\n', ...
         i, qq(i,X), qq(i,Y), qq(i,TH));
   end
end

%% Dibujar la estructura y su deformada
esc_def    = 50;            % escalamiento de la deformada
esc_faxial = 0.2;           % escalamiento del diagrama de axiales
esc_V      = 0.3;           % escalamiento del diagrama de cortantes
esc_M      = 0.3;           % escalamiento del diagrama de momentos

xdef = xnod + esc_def*vect_mov(:,[X Y]);

figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
figure(3); hold on; title('Fuerza axial [ton]');      xlabel('x, m'); ylabel('y, m'); axis equal
figure(4); hold on; title('Fuerza cortante [ton]');   xlabel('x, m'); ylabel('y, m'); axis equal
figure(5); hold on; title('Momento flector [ton-m]'); xlabel('x, m'); ylabel('y, m'); axis equal

for e = 1:nbar
   x1 = xnod(LaG(e,NL1), X);  x2 = xnod(LaG(e,NL2), X);
   y1 = xnod(LaG(e,NL1), Y);  y2 = xnod(LaG(e,NL2), Y);

   dibujar_barra_deformada_portico(...
      A(mat(e)), E(mat(e)), I(mat(e)), ...
      x1,y1, x2,y2, qxloc{e}, qyloc{e}, ...
      qe_loc{e}, T{e}*a(idx(e,:)), ...
      esc_def, esc_faxial, esc_V, esc_M);
end

%% bye, bye!
