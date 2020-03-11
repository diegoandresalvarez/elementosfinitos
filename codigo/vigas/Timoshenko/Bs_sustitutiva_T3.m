% Imposicion de las deformaciones angulares gxz
% Se realiza para el elemento finito de Timoshenko de tres nodos

syms G Aast L xi w1 w2 w3 t1 t2 t3

% La matriz Bs "normal"
Bs = @(xi) (2/L)*[ xi - 1/2, -(L*xi*(xi - 1))/4, -2*xi, ...
                           (L*(xi^2 - 1))/2, xi + 1/2, -(L*xi*(xi + 1))/4];

% Las funciones de interpolacion de las deformaciones gxz                        
Ng1 = poly2sym(polyfit([-1/sqrt(3) +1/sqrt(3)], [1 0], 1), xi); 
Ng2 = poly2sym(polyfit([-1/sqrt(3) +1/sqrt(3)], [0 1], 1), xi); 

% Se calcula la matriz Bc_sustitutiva
Bs_sustitutiva = Ng1*Bs(sym(-1/sqrt(3))) + Ng2*Bs(sym(+1/sqrt(3)));

% y la correspondiente matriz de rigidez
Ks = int(Bs_sustitutiva.'*G*Aast*Bs_sustitutiva*L/2,xi,-1,1);

disp('Ks = ((G*Aast)/(9*L)) * '), disp(Ks/(G*Aast/(9*L)))

% Este Ks coincide con aquel integrado con dos puntos de GL.

% bye bye!
