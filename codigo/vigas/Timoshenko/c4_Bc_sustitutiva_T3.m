% Imposicion de las deformaciones angulares gxz
% Se realiza para el elemento finito de Timoshenko de tres nodos

syms G Aast L xi w1 w2 w3 t1 t2 t3

% La matriz Bc "normal"
Bc = @(xi) (2/L)*[ xi - 1/2, -(L*xi*(xi - 1))/4, -2*xi, ...
                           (L*(xi^2 - 1))/2, xi + 1/2, -(L*xi*(xi + 1))/4];

% Las funciones de interpolacion de las deformaciones gxz                        
Ng1 = poly2sym(polyfit([-1/sqrt(3) +1/sqrt(3)], [1 0], 1), xi); 
Ng2 = poly2sym(polyfit([-1/sqrt(3) +1/sqrt(3)], [0 1], 1), xi); 

% Se calcula la matriz Bc_sustitutiva
Bc_sustitutiva = Ng1*Bc(sym(-1/sqrt(3))) + Ng2*Bc(sym(+1/sqrt(3)));

% y la correspondiente matriz de rigidez
Kc = int(Bc_sustitutiva.'*G*Aast*Bc_sustitutiva*L/2,xi,-1,1);

disp('Kc = ((G*Aast)/(9*L)) * '), disp(Kc/(G*Aast/(9*L)))

% Este Kc coincide con aquel integrado con dos puntos de GL.

