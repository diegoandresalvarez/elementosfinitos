clear, clc

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta 
gs = sym('gs', [1 6]).';

%% Se definen los puntos de colocacion
d = sym(sqrt(3))/6;

nod = [ ...
% xi      eta       beta
 0.5-d     0           0    % 1
 0.5+d     0           0    % 2
 0.5+d   0.5-d    3*pi/4    % 3
 0.5-d   0.5+d    3*pi/4    % 4
   0     0.5+d      pi/2    % 5
   0     0.5-d      pi/2 ]; % 6

n = size(nod,1);           % numero de puntos de colocacion
beta = nod(:,3);           % se definen las direcciones

%% Se calculan los componentes de gs
gxi = gs.*cos(beta);        geta  = gs.*sin(beta);
gpg = [ gxi; geta ]; 

%% Se calculan las ecuaciones que forman el sistema de ecuaciones 9.57
xxi  = nod(:, XI); 
eeta = nod(:, ETA);

A = blkdiag([ ones(n,1) xxi eeta ],[ ones(n,1) xxi eeta ]);
a = sym('a',[1 6]).';
igualdad = A*a == gpg

IG = [
igualdad([1 2])
igualdad(3)*cos(beta(3)) + igualdad( 9)*sin(beta(3))  %se forza sin(t^2) + cos(t^2) = 1
igualdad(4)*cos(beta(4)) + igualdad(10)*sin(beta(4))
igualdad([11 12])
]

%% Se calcula la matriz P
syms a1 a2 a3 a4 a5 a6
ee = mat2cell(eye(6),[ 1 1 1 1 1 1 ], 6 ); % produce {[1 0 0 0 0 0];[0 1 0 0 0 0 ]; ... ;[0 0 0 0 0 1]}
P = subs(lhs(IG), {a1, a2, a3, a4, a5, a6}, ee)

%% Se calcula la matriz T
TT = cell(1,n);
for i = 1:n
   TT{i} = [ cos(beta(i)) sin(beta(i)) ];
end
T = blkdiag(TT{:})                              % eq 6.65, 6.89

%% Se define la matriz A
A = [ 1  xi eta   0   0   0 
      0  0   0   1  xi eta ];   % eq 6.62, 6.116

A_invP_T = A*inv(P)*T

invP_T_metodo1 = double(P\T)
