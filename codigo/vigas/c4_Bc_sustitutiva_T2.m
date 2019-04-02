% Imposicion de las deformaciones angulares gxz
% Se realiza para el elemento finito de Timoshenko de dos nodos

syms G Aast L xi w1 w2 t1 t2

% La matriz Bc "normal"
Bc = @(xi) [ -1/L, xi/2 - 1/2, 1/L, - xi/2 - 1/2];

% Las funciones de interpolacion de las deformaciones gxz                        
Ng = 1; % se evaluara en el centro

% Se calcula la matriz Bc_sustitutiva
Bc_sustitutiva = Ng*Bc(0);

% y la correspondiente matriz de rigidez
Kc = int(Bc_sustitutiva.'*G*Aast*Bc_sustitutiva*L/2,xi,-1,1);

disp('Kc = (G*Aast/L) * '), disp(Kc/(G*Aast/L))

% Este Kc coincide con aquel integrado con un punto de GL.
