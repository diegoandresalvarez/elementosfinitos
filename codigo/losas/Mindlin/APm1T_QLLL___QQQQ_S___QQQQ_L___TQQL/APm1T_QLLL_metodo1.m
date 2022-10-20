clear, clc, close all

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta

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

% se crean las variables gpg
ngamma = size(idx,2);  % numero de puntos de colocación por dirección
gpg = [ crear_gpg('xi',  idx(XI, :))    % Los + = gxi1  gxi3
        crear_gpg('eta', idx(ETA,:)) ]; % Los x = geta2 geta4

% se calculan las funciones de forma en cada punto de colocación
N = sym(zeros(2,ngamma));
for i = [ XI ETA ]
   xxi  = nod(idx(i,:), XI);
   eeta = nod(idx(i,:), ETA);
   switch i
      case XI
         A         = [ ones(ngamma,1) eeta ];
         variables = [ 1 eta ];
      case ETA
         A         = [ ones(ngamma,1) xxi ];
         variables = [ 1 xi ];
   end

   for j = 1:ngamma
      % se arma el sistema de ecuaciones
      b = zeros(ngamma,1);   b(j) = 1;
      coef_alpha = A\b;
      N(i,j) = simplify(variables*coef_alpha);
      fprintf('j = %d (%s): %s\n', j, char(gpg(i,j)), char(N(i,j)));
   end
   fprintf('----------------------------------------------------------\n');
end

%% Se imprime el polinomio de interpolacion
gp = sum(simplify(N.*gpg), 2)

%% Se calcula la matriz A*inv(P)*T
npc_max = max(max(idx));
gpg = [ crear_gpg('xi',  1:npc_max)    % Los + = gxi1  .. gxi4
        crear_gpg('eta', 1:npc_max) ]; % Los x = geta1 .. geta4
gpg = gpg(:).';  % gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4

A_invP_T = sym(zeros(2, 2*npc_max));
for i = 1:2*npc_max
   A_invP_T(XI, i) = feval(symengine, 'coeff', gp(1), gpg(i), 1);
   A_invP_T(ETA,i) = feval(symengine, 'coeff', gp(2), gpg(i), 1);   
end
gp_metodo1 = simplify(A_invP_T*gpg.')

A_invP_T_metodo1 = A_invP_T

%% Graficamos las funciones de forma en M=A_invP_T
graficar_A_invP_T(nod, idx, A_invP_T)
