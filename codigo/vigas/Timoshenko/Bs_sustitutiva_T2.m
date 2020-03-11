% Imposicion de las deformaciones angulares gxz
% Se realiza para el elemento finito de Timoshenko de dos nodos

syms G Aast L xi w1 w2 t1 t2

% La matriz Bs "normal"
Bs = @(xi) [ -1/L, xi/2 - 1/2, 1/L, - xi/2 - 1/2];

% Las funciones de interpolacion de las deformaciones gxz                        
Ng = 1; % se evaluara en el centro

% Se calcula la matriz Bs_sustitutiva
Bs_sustitutiva = Ng*Bs(0);

% y la correspondiente matriz de rigidez
Ks = int(Bs_sustitutiva.'*G*Aast*Bs_sustitutiva*L/2,xi,-1,1);

disp('Ks = (G*Aast/L) * '), disp(Ks/(G*Aast/L))

% Este Ks coincide con aquel integrado con un punto de GL.

% bye bye!
