clear,clc        % borra la memoria y la pantalla

%% Programa para deducir la matriz de rigidez para un elemento rectangular 
% de 4 nodos

%% Se definen las variables simbolicas
syms a b r s d11 d12 d21 d22 d33 t a1 a2 a3 a4 a5 a6

d21 = d12;
d11 = d22;

%% Funciones de forma del rectangulo
N1 = (1-r/a)*(1-s/b)/4;
N2 = (1+r/a)*(1-s/b)/4;
N3 = (1+r/a)*(1+s/b)/4;
N4 = (1-r/a)*(1+s/b)/4;

%% se ensambla la matriz de funciones de forma
N = [ N1 0  N2 0  N3 0  N4 0
      0  N1 0  N2 0  N3 0  N4 ];

%% matriz constitutiva (ya luego uno particularizara si es para tension o
% para deformacion plana)
D = [d11 d12 0
     d21 d22 0
     0   0   d33];
  
%% matriz de deformaciones
B1 = [diff(N1,r)   0         
      0            diff(N1,s)
      diff(N1,s)   diff(N1,r)];
      
B2 = [diff(N2,r)   0
      0            diff(N2,s)
      diff(N2,s)   diff(N2,r)];
      
B3 = [diff(N3,r)   0
      0            diff(N3,s)
      diff(N3,s)   diff(N3,r)];
      
B4 = [diff(N4,r)   0
      0            diff(N4,s)
      diff(N4,s)   diff(N4,r)];

% se ensambla la matriz de deformaciones
B = [B1 B2 B3 B4];

%% realizo la integracion
K = int(int(B.'*D*B*t, r, -a,a), s,-b,b);
K = expand(K);

% Sustitucion de variables
K = subs(K,t*b*d11/(6*a),a1);    K = subs(K,-t*b*d11/(6*a),-a1);
K = subs(K,t*a*d22/(6*b),a2);    K = subs(K,-t*a*d22/(6*b),-a2);
K = subs(K,d12*t/4,      a3);    K = subs(K,-d12*t/4,      -a3);
K = subs(K,a*d33*t/(6*b),a4);    K = subs(K,-a*d33*t/(6*b),-a4);
K = subs(K,t*b*d33/(6*a),a5);    K = subs(K,-t*b*d33/(6*a),-a5);
K = subs(K,t*d33/4,      a6);    K = subs(K,-t*d33/4,      -a6);

K = subs(K,t*b*d11/(3*a),2*a1);  K = subs(K,-t*b*d11/(3*a),-2*a1);
K = subs(K,t*a*d22/(3*b),2*a2);  K = subs(K,-t*a*d22/(3*b),-2*a2);
K = subs(K,d12*t/2,      2*a3);  K = subs(K,-d12*t/2,      -2*a3);
K = subs(K,a*d33*t/(3*b),2*a4);  K = subs(K,-a*d33*t/(3*b),-2*a4);
K = subs(K,t*b*d33/(3*a),2*a5);  K = subs(K,-t*b*d33/(3*a),-2*a5);
K = subs(K,t*d33/2,      2*a6);  K = subs(K,-t*d33/2,      -2*a6);

syms a14 a25 a36 c41 c14 b36 b63 c25 c52
K = subs(K,a1 + a4, a14);        K = subs(K,-a1 - a4, -a14);
K = subs(K,a2 + a5, a25);        K = subs(K,-a2 - a5, -a25);
K = subs(K,a2 - 2*a5, c25);      K = subs(K,a5 - 2*a2, c52);
K = subs(K,a4 - 2*a1, c41);      K = subs(K,a1 - 2*a4, c14);
K = subs(K,a3 - a6,   b36);      K = subs(K,a6 - a3,   b63);
K = subs(K,a3 + a6, a36);        K = subs(K,-a3 - a6, -a36);

%% finalmente se imprime la matriz
K

disp(' ');
disp('Para facilitar la comparacion de la matriz con el libro de Onate')
disp('se muestra unicamente la parte superior de K =')
triu(K)

%%
%% se calcula el vector de fuerzas nodales equivalentes para las fuerzas masicas
syms bx by
bb = [bx; by];

fb = int(int(N.'*bb*t, r, -a,a), s,-b,b);
disp('K = a*b*t*');
pretty(fb/(a*b*t));

%% Se calcula la matriz de masa consistente
syms rho
M = t*rho*int(int(N.'*N, r, -a,a), s,-b,b)
disp('M = rho*t*Area/36 * ')
disp(36*M/(rho*t*4*a*b))