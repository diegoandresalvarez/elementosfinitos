% Se definen las variables simbolicas
syms b E A L u3 u4 P

% Se calculan las rigideces para cada barra
k1 = E*A/L; k2 = E*A/L; k3 = 2*E*A/L;

% Se define la matriz de rigidez
K = [ ...
     k1   0   -k1         0
     0    k2  -k2         0
    -k1  -k2   k1+k2+k3  -k3
     0    0   -k3         k3 ]

% Se define el vector de desplazamientos nodales, teniendo en cuenta que 
% u1=0 y u2=0
a = [0; 0; u3; u4]

% Se define el vector de fuerzas nodales de equilibrio
%f = [ 0; 0; 0; P]
f = [ b*L/2; b*L/2; P/2 + b*L; P]

% Se definen g.d.l. conocidos y desconocidos asociados a los desplazamientos
c = [1 2];         d = [3 4];

% Se descomponen los vectores a, f y la matriz K
Kcc = K(c,c);      Kcd = K(c,d);
Kdc = K(d,c);      Kdd = K(d,d);
ac  = a(c);      % ad  = a(d);
fc  = f(d);        fd  = f(c);

% Se calculan los vectores ad y fd
ad = Kdd\(fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad - fd

