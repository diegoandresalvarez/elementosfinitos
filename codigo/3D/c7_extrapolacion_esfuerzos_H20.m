clear, clc
[x_gl, w_gl] = gausslegendre_quad_hexa(2);
n_gl = length(w_gl);

A1 = zeros(n_gl);
for i = 1:n_gl
   xi_gl   = x_gl(i,1);
   eta_gl  = x_gl(i,2);
   zeta_gl = x_gl(i,3);
   A1(i,:) = [ 1 xi_gl eta_gl zeta_gl xi_gl*eta_gl xi_gl*zeta_gl eta_gl*zeta_gl xi_gl*eta_gl*zeta_gl];
end

nod = [ ...
%  xi   eta  zeta   nodo
   -1   -1   -1   %  1
    0   -1   -1   %  2
    1   -1   -1   %  3
    1    0   -1   %  4
    1    1   -1   %  5
    0    1   -1   %  6
   -1    1   -1   %  7
   -1    0   -1   %  8
   -1   -1    0   %  9
    1   -1    0   % 10
    1    1    0   % 11
   -1    1    0   % 12
   -1   -1    1   % 13
    0   -1    1   % 14
    1   -1    1   % 15
    1    0    1   % 16
    1    1    1   % 17
    0    1    1   % 18
   -1    1    1   % 19
   -1    0    1 ];% 20

A2 = zeros(20,8);
for i = 1:20
      xi_gl   = nod(i,1);
      eta_gl  = nod(i,2);
      zeta_gl = nod(i,3);
      A2(i,:) = [ 1 xi_gl eta_gl zeta_gl xi_gl*eta_gl xi_gl*zeta_gl eta_gl*zeta_gl xi_gl*eta_gl*zeta_gl];
end

A1 = sym(A1);
A2 = sym(A2);

A = A2/A1 %= A2*inv(A1)
