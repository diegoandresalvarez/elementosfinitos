clear, clc, close all

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta gxi1 gxi2 gxi3 gxi4 gxi5 geta3 geta6 geta7 geta8 geta9
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
     
gpg = [ gxi1  gxi2  gxi3  gxi4  gxi5      % Los +
        geta3 geta6 geta7 geta8 geta9 ];  % Los x   
     
ngamma = size(idx,2);  % numero de puntos que definen las funciones de forma
N = sym(zeros(2,ngamma));

for i = [ XI ETA ]
   xxi  = nod(idx(i,:), XI); 
   eeta = nod(idx(i,:), ETA);
   switch i
      case XI
         A         = [ ones(ngamma,1) xxi eeta  xxi.*eeta  eeta.^2 ];
         variables = [ 1 xi eta xi*eta eta^2 ];
      case ETA
         A         = [ ones(ngamma,1) xxi eeta  xxi.*eeta  xxi.^2 ];      
         variables = [ 1 xi eta xi*eta xi^2 ];
   end

   for j = 1:ngamma
      % se arma el sistema de ecuaciones
      b = zeros(ngamma,1);   b(j) = 1;
      coef_alpha = A\b;
      fprintf('j = %d (%s):  ', j, char(gpg(i,j))); 
      N(i,j) = simplify(variables*coef_alpha);
      disp(N(i,j))
   end
   fprintf('----------------------------------------------------------\n');
end

%% Se imprime el polinomio de interpolacion
gp = sum(simplify(N.*gpg), 2)

syms gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4
syms gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9

%% Se calcula la matriz A*inv(P)*T
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 ...
        gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9 ];

A_invP_T = sym(zeros(2, 18));

for i = 1:length(gpg)
   A_invP_T(XI, i) = feval(symengine, 'coeff', gp(1), gpg(i), 1);
   A_invP_T(ETA,i) = feval(symengine, 'coeff', gp(2), gpg(i), 1);   
end
gp_metodo1 = simplify(A_invP_T*gpg.')

A_invP_T_metodo1 = A_invP_T

%% Graficamos las funciones de forma en M=A_invP_T
graficar_A_invP_T(nod, idx, A_invP_T)