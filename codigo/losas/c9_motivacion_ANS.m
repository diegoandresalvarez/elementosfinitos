% Funciones de forma del rectangulo de 4 nodos
syms x y a b

% Se definen las funciones de forma
N = sym(zeros(4,1));
N(1) = (1 - x/a)*(1 - y/b)/4;
N(2) = (1 + x/a)*(1 - y/b)/4;
N(3) = (1 + x/a)*(1 + y/b)/4;
N(4) = (1 - x/a)*(1 + y/b)/4;

% Se definen los vectores de desplazamientos y giros nodales
we  = sym('w',  [4 1]);   w  = sum(N.*we);
txe = sym('tx', [4 1]);   tx = sum(N.*txe);
tye = sym('ty', [4 1]);   ty = sum(N.*tye);

% Se calculan las deformaciones gxz y se extraen sus coeficientes
gxz = diff(w,x) - tx;

a = sym(zeros(4,1));
a(1) = feval(symengine, 'coeff', gxz, '[x,y]', '[0,0]');
a(2) = feval(symengine, 'coeff', gxz, '[x,y]', '[0,1]');
a(3) = feval(symengine, 'coeff', gxz, '[x,y]', '[1,0]');
a(4) = feval(symengine, 'coeff', gxz, '[x,y]', '[1,1]');

condicion = simple(a == 0)
condicion(3) - condicion(4)  % restando las condiciones 3 y 4

