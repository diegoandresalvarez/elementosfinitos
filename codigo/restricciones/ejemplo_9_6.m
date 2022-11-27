clear, clc

% Onate, Vol 1, pag 317.
% Calcule los desplazamientos en los nodos asumiendo que u1 == u2
%
% Modulo de elasticidad E
% Area transversal A
%                                         |/
%   P=>         P=>         P=>           |/
%     1-----------2-----------3-----------4/
%                                         |/
%     |<--- L --->|<--- L --->|<--- L --->|/
%

% Se definen las variables simbolicas
syms u1 u2 u3 P alpha k

% Se calculan las rigideces para cada barra
k1 = k; k2 = k; k3 = k;

% Se define la matriz de rigidez
K = [ k1  -k1       0        0
     -k1   k1+k2   -k2       0
      0   -k2       k2+k3   -k3
      0    0       -k3       k3 ];

% Se define el vector de desplazamientos nodales, teniendo en cuenta que 
% u4=0
a = [u1; u2; u3; 0];

% Se define el vector de fuerzas nodales de equilibrio
f = [ P; P; P; 0 ];

% IMPONDREMOS LA RESTRICCION u1 == u2 y se coloca en C*a = g
r = [ 1 3 4 ];     % nodos a retener (incluye el nodo maestro = 1)
e = [ 2 ];         % nodo esclavo
c = [ 4 ];         % GDL del empotramiento
C = [ 1 -1 0 0 ];
g = 0;

n       = length(a);  % numero de GDL de la estructura
num_res = size(C,1);  % numero de GDL esclavos o restringidos = length(e)
n_c = n - num_res;    % numero de GDL no restringidos         = length(r)

%% METODO 1: transformacion de K*a-f = q
%{
Cr = C(:,r);   Ce = C(:,e);

% Se verifica si Ce es singular o no
if abs(det(Ce)) < 1e-5
    error('Ce debe ser invertible');
end

I  = eye(n_c);
T  = [  I
       -Ce\Cr ]; % = -inv(Ce)*Cr

O  = zeros(n_c, 1);
g0 = [  O
        Ce\g  ]; % =  inv(Ce)*g

% se reordenan los GDL de K y de f
idx_re = [r e];

Kre = K(idx_re, idx_re);
fre = f(idx_re);
   
Kr = T'*Kre*T;
fr = T'*(fre - Kre*g0);

% Se definen g.d.l. conocidos y desconocidos asociados a los desplazamientos
% idx con la numeracion original (del grafico)
% obseve que en Kr*ar - fr = qr no hay nodos esclavos
d = setdiff(r, c);

% idx con la numeracion del reordenamiento [r e], esto es idx_re
c_re = arrayfun(@(x) find(r==x), c);
d_re = arrayfun(@(x) find(r==x), d);

% Se descomponen los vectores a, f y la matriz K
Krcc = Kr(c_re,c_re);      Krcd = Kr(c_re,d_re);
Krdc = Kr(d_re,c_re);      Krdd = Kr(d_re,d_re);

ar   = a(r);
arc  = ar(c_re);

frc  = fr(d_re);           frd  = fr(c_re);

% Se calculan los vectores ard y qrd (de los GDL a retener)
ard = Krdd\(frc - Krdc*arc);
qrd = Krcc*arc + Krcd*ard - frd;

% se calcula el vector a_re
ar = sym(zeros(n_c,1));
ar(c_re) = arc;
ar(d_re) = ard;
a_re = T*ar + g0;

% se calcula el vector q_re
q_re = sym(zeros(n,1));
q_re(c_re) = qrd;

% se reordena a de la indexacion re a la indexacion original (del grafico)
idx_orig = arrayfun(@(x) find(idx_re==x), 1:n);
a = a_re(idx_orig);
q = q_re(idx_orig);

disp('a = '); pretty(a);
disp('q = '); pretty(q);
%}

%% METODO 2: multiplicadores de Lagrange
O  = zeros(num_res);
Kr = [ K C'
       C O ];
fr = [ f
       g ];
   
% se definen los GDL conocidos y desconocidos asociados a los desplazamientos
d = setdiff(1:(n + num_res), c);

lambda = sym(zeros(num_res,1));
ar = [ a
       lambda ];

% se descomponen los vectores a, f y la matriz K
Krcc = Kr(c,c);      Krcd = Kr(c,d);
Krdc = Kr(d,c);      Krdd = Kr(d,d);
arc  = ar(c);
frc  = fr(d);        frd  = fr(c);

% se calculan los vectores ard y frd (de los GDL a retener)
ard = Krdd\(frc - Krdc*arc);
qrd = Krcc*arc + Krcd*ard - frd;

% se calcula de nuevo el vector a
ar     = sym(zeros(n + num_res, 1));
ar(c)  = arc;
ar(d)  = ard;
a      = ar(1:n);

% Observe que lambda dara una fuerza, por lo que los multiplicadores de 
% Lagrange se pueden interpretar como las fuerzas de enlace necesarias para
% imponer las restricciones en los GDL.
lambda = ar(n+1:n+num_res)

fuerzas_de_enlace = -C'*lambda


% se calcula de nuevo el vector q
qr     = sym(zeros(n + num_res, 1));
qr(c)  = qrd;
q      = qr(1:n);

disp('a = '); pretty(a);
disp('q = '); pretty(q);
%}

%% METODO 3: metodo de las penalizaciones
%{
% se definen los GDL conocidos y desconocidos asociados a los desplazamientos
d = setdiff(1:n, c);

% alpha = diag([ a1, a2, ..., ar ]);
% Observe que C'*alpha*C se puede interpretar como la matriz de rigidez
asociada a una barra con rigidez alpha que une los nodos 1 y 2
Kb = K + C'*alpha*C;
fb = f + C'*alpha*g;

% Se descomponen los vectores a, f y la matriz K
Kbcc = Kb(c,c);      Kbcd = Kb(c,d);
Kbdc = Kb(d,c);      Kbdd = Kb(d,d);
ac   = a(c); 
fbc  = fb(d);        fbd  = fb(c);

% Se calculan los vectores ad y fd
ad = Kbdd\(fbc - Kbdc*ac);
qd = Kbcc*ac + Kbcd*ad - fbd;

% se calcula de nuevo el vector a
a    = sym(zeros(n, 1));
a(c) = ac;
a(d) = ad;

% se calcula de nuevo el vector q
q    = sym(zeros(n, 1));
q(c) = qd;

pretty(partfrac(a, k)) % presentemos el resultado como fracciones parciales
disp('a = '); pretty(limit(a, alpha, inf))
disp('q = '); pretty(q);
%}

%% Bye, bye!

