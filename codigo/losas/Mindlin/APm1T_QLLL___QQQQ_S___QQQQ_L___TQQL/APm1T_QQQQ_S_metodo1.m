clear, clc, close all

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta
a = 1/sqrt(sym(3));

%% Se definen los puntos de colocacion
nod = [ ...
% xi eta
   a   1    % 1
  -a   1    % 2
   0   0    % 3
   a  -1    % 4
  -a  -1    % 5
   1   a    % 6
   1  -a    % 7
  -1   a    % 8
  -1  -a    % 9
];

%% Se calculan las funciones de forma bidimensionales
idx = [ 1 2 3 4 5    % puntos usados por gamma_xi  +
        3 6 7 8 9 ]; % puntos usados por gamma_eta x

% se crean las variables gpg
ngamma_1_dir = size(idx,2);  % numero de puntos de colocación por dirección
gpg = [ crear_gpg('xi',  idx(XI, :))    % Los + = gxi1  gxi2  gxi3  gxi4  gxi5
        crear_gpg('eta', idx(ETA,:)) ]; % Los x = geta3 geta6 geta7 geta8 geta9

% se calculan las funciones de forma en cada punto de colocación
N = sym(zeros(2,ngamma_1_dir));
for i = [ XI ETA ]
   xxi  = nod(idx(i,:), XI);
   eeta = nod(idx(i,:), ETA);
   switch i
      case XI
         A         = [ ones(ngamma_1_dir,1) xxi eeta  xxi.*eeta  eeta.^2 ];
         variables = [ 1 xi eta xi*eta eta^2 ];
      case ETA
         A         = [ ones(ngamma_1_dir,1) xxi eeta  xxi.*eeta  xxi.^2 ];
         variables = [ 1 xi eta xi*eta xi^2 ];
   end

   for j = 1:ngamma_1_dir
      % se arma el sistema de ecuaciones
      b = zeros(ngamma_1_dir,1);   b(j) = 1;
      coef_alpha = A\b;
      N(i,j) = simplify(variables*coef_alpha);
      fprintf('j = %d (%s): %s\n', j, char(gpg(i,j)), char(N(i,j)));
   end
   fprintf('----------------------------------------------------------\n');
end

%% Se imprime el polinomio de interpolacion
gp = sum(simplify(N.*gpg), 2)

%% Se calcula la matriz A*inv(P)*T
ngamma = max(max(idx));
gpg = [ crear_gpg('xi',  1:ngamma)    % Los + = gxi1  .. gxi9
        crear_gpg('eta', 1:ngamma) ]; % Los x = geta1 .. geta9
gpg = gpg(:).';  % gxi1 geta1 gxi2 geta2 .. gxi9 geta9

A_invP_T = sym(zeros(2, 2*ngamma));
for i = 1:2*ngamma
   A_invP_T(XI, i) = feval(symengine, 'coeff', gp(1), gpg(i), 1);
   A_invP_T(ETA,i) = feval(symengine, 'coeff', gp(2), gpg(i), 1);   
end
gp_metodo1 = simplify(A_invP_T*gpg.')

A_invP_T_metodo1 = A_invP_T

%% Graficamos las funciones de forma en M=A_invP_T
graficar_A_invP_T(nod, idx, A_invP_T)
