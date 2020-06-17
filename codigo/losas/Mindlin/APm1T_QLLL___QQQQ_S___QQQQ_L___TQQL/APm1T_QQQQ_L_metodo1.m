clear
clc

%% Se definen las constantes
XI = 1; ETA = 2;
syms xi eta gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 gxi6 geta6 
syms        gxi7 geta7 gxi8 geta8 gxi9 geta9 gxi10 geta10 gxi11 geta11 gxi12 geta12

a = 1/sym(sqrt(3));

%% Se definen los puntos de colocacion
nod = [ ...
% xi eta
   a   1    %  1
  -a   1    %  2
   a   0    %  3
  -a   0    %  4
   a  -1    %  5 
  -a  -1    %  6
   1   a    %  7
   1  -a    %  8          
   0   a    %  9
   0  -a    % 10          
  -1   a    % 11
  -1  -a ]; % 12

%% Se calculan las funciones de forma bidimensionales
idx = [ 1 2 3 4 5 6      % puntos usados por gamma_xi  +
        7 8 9 10 11 12]; % puntos usados por gamma_eta x
     
gpg = [ gxi1  gxi2  gxi3  gxi4   gxi5   gxi6      % Los +
        geta7 geta8 geta9 geta10 geta11 geta12];  % Los x   
     
n = size(idx,2);  % numero de puntos que definen las funciones de forma
N = sym(zeros(2,n));

for i = 1:2
   xxi  = nod(idx(i,:), XI); 
   eeta = nod(idx(i,:), ETA);
   switch i
      case 1
         A         = [ ones(n,1) xxi eeta  xxi.*eeta  eeta.^2  xxi.*eeta.^2];
         variables = [ 1  xi  eta  xi*eta  eta^2  xi*eta^2 ];
      case 2
         A         = [ ones(n,1) xxi eeta  xxi.*eeta  xxi.^2 eeta.*xxi.^2];      
         variables = [ 1 xi eta xi*eta xi^2 eta*xi^2];
   end

   for j = 1:n
      % se arma el sistema de ecuaciones
      b = zeros(n,1);   b(j) = 1;
      coef_alpha = A\b;
      fprintf('j = %d (%s):  ', j, char(gpg(i,j))); 
      N(i,j) = simplify(variables*coef_alpha);
      disp(N(i,j))
   end
   fprintf('----------------------------------------------------------\n');
end

%% Se imprime el polinomio de interpolacion
gp = sum(simplify(N.*gpg), 2)

%% Se calcula la matriz A*inv(P)*T
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 ...
        gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9 ...
        gxi10 geta10 gxi11 geta11 gxi12 geta12 ];

A_invP_T = sym(zeros(2, 24));

for i = 1:length(gpg)
   A_invP_T(1,i) = feval(symengine, 'coeff', gp(1), gpg(i), 1);
   A_invP_T(2,i) = feval(symengine, 'coeff', gp(2), gpg(i), 1);   
end
gp_metodo1 = simplify(A_invP_T*gpg.')

A_invP_T_metodo1 = A_invP_T