%% Deduccion de las funciones de forma del elemento de barra de 2 nodos

syms u1 u2 x x1 x2 a1 a0;
r = solve(u1 == a1*x1 + a0, ...
          u2 == a1*x2 + a0,     a0, a1);
disp('a0 = '); pretty(r.a0)
disp('a1 = '); pretty(r.a1)

u = r.a1*x + r.a0;        % Se define ahora u(x) ya que conocemos a1 y a0
u = collect(u, [u1, u2]); % Se factoriza u1 y y2
disp('u = '); pretty(u)   % Observe aqui las funciones de forma