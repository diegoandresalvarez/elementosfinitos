clear, clc, close all

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4

%          ^ eta
%          |
%          |
%          3
%   (4)----+----(3)
%    |           |  
%    |           |
%  4 x           x 2  ------> xi
%    |           |
%    |           |
%   (1)----+----(2)
%          1

%% Se definen los puntos de colocacion
nod = [...
% xi  eta
   0   -1   % 1
   1    0   % 2
   0    1   % 3
  -1    0]; % 4

%% Se calculan las funciones de forma bidimensionales
idx = [ 1 3    % puntos usados por gamma_xi  +
        2 4 ]; % puntos usados por gamma_eta x
     
gpg = [ gxi1  gxi3      % Los +
        geta2 geta4 ];  % Los x   
     
n = size(idx,2);  % numero de puntos que definen las funciones de forma
N = sym(zeros(2,n));

for i = 1:n
   xxi  = nod(idx(i,:), XI); 
   eeta = nod(idx(i,:), ETA);
   A         = [ ones(n,1) xxi eeta ];
   variables = [ 1         xi  eta  ];

   for j = 1:n
      % se arma el sistema de ecuaciones
      b = zeros(n,1);   b(j) = 1;
      coef_alpha = A\b;
      %fprintf('j = %d (%s):  ', j, char(gpg(i,j))); 
      N(i,j) = simplify(variables*coef_alpha);
      %disp(N(i,j))
   end
   %fprintf('----------------------------------------------------------\n');
end

%% Se imprime el polinomio de interpolacion
gp = sum(simplify(N.*gpg), 2)

%% Se calcula la matriz A*inv(P)*T
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 ];

A_invP_T = sym(zeros(2, 8)); % esta es la matriz M (paso 2, metodo 1)

for i = 1:length(gpg)
   A_invP_T(1,i) = feval(symengine, 'coeff', gp(1), gpg(i), 1);
   A_invP_T(2,i) = feval(symengine, 'coeff', gp(2), gpg(i), 1);   
end
gp_metodo1 = simplify(A_invP_T*gpg.')

A_invP_T_metodo1 = A_invP_T

%% Graficamos las funciones de forma en M=A_invP_T
graficar_A_invP_T(nod, idx, A_invP_T)