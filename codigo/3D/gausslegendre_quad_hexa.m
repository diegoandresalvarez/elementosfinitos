function [x_gl, w_gl] = gausslegendre_quad_hexa(n_gl)

[x, w] = gausslegendre_quad(n_gl);

x_gl = zeros(n_gl^2,3);
w_gl = zeros(n_gl^2,1);

r = 0;
for i = 1:n_gl
   for j = 1:n_gl
      for k = 1:n_gl
         r = r+1;
         x_gl(r,:) = [ x(i) x(j) x(k) ];
         w_gl(r) = w(i)*w(j)*w(k);
      end
   end
end

return;