%% Interpolacion acoplada

syms xi w1 w2 w3 t1 t2 t3 L

dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Se definen w y t (w3, t3 es el nodo central)
% la idea de esta numeracion es que eliminaremos el nodo central y lo
% expresaremos en funcion de los gdl asociados a los nodos de los lados
w = N1*w1 + N2*w3 + N3*w2;  % ¡ojo al orden de los terminos!
t = N1*t1 + N2*t3 + N3*t2;  % ¡ojo al orden de los terminos!

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = simplify(diff(w,xi)*dxi_dx - t)
A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);
% Nota: tambien se pudo haber usado:
%{
gxz = simplify(diff(w,xi)*dxi_dx - t)
tmp = coeffs(gxz, xi);
A = tmp(1);
B = tmp(2);
C = tmp(3);
%}

%% Los terminos que acompanan a xi^2 y xi deben hacerse igual a cero (B = C = 0)
sol = solve(B==0, C==0, w3,t3);
disp('w3 = '); disp(sol.w3)
disp('t3 = '); disp(sol.t3)

%% Y en funcion de w3 y t3 se reescriben w y t
w = N1*w1 + N2*sol.w3 + N3*w2;
w = collect(w, {'w1', 'w2'})

t = N1*t1 + N2*sol.t3 + N3*t2;
t = collect(t, {'t1', 't2'})

%% Se recalcula gxz
gxz = simplify(diff(w,xi)*dxi_dx - t)

%% bye, bye!
return
