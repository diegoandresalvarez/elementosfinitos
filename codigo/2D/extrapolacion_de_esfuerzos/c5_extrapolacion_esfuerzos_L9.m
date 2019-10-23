n_gl = 2;
[x_gl, w_gl]  = gausslegendre_quad(n_gl);
A1 = zeros(n_gl^2);

i = 0;
%pos = [1 2 4 3];
pos = [1 2 3 4];  % programa
for p = 1:n_gl
   for q = 1:n_gl
      i = i+1;
      xi_gl  = x_gl(p);
      eta_gl = x_gl(q);  
      A1(pos(i),:) = [ 1 xi_gl eta_gl xi_gl*eta_gl ];
   end
end

% Coordenadas de los nodos
%
% Numeracion local:
%      ^ eta
%      |
%      |
%  7---6---5
%  |   |   |
%  8---9---4----> xi
%  |   |   |
%  1---2---3

nod = [ ...
%  xi   eta     % nodo   
   -1   -1      %  1
    0   -1      %  2
    1   -1      %  3
    1    0      %  4
    1    1      %  5
    0    1      %  6
   -1    1      %  7
   -1    0      %  8
    0    0  ];  %  9

A2 = zeros(9,4);
for i = 1:9
      xi_gl  = nod(i,1);
      eta_gl = nod(i,2);
      
      A2(i,:) = [ 1 xi_gl eta_gl xi_gl*eta_gl ];      
end

A1 = sym(A1);
A2 = sym(A2);

A = A2/A1 %= A2*inv(A1)
