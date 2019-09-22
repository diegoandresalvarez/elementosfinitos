function z = int_f(N,x1,x2,x3,y1,y2,y3)

%   This function evaluates \iint_K f(x,y) dxdy using 
%   the Gaussian quadrature of order N where K is a 
%   triangle with vertices (x1,y1), (x2,y2) and (x3,y3).

% Tomada de: http://math2.uncc.edu/~shaodeng/TEACHING/math5172/2010Spring/

xw = TriGaussPoints(N);  % get quadrature points and weights 

% calculate the area of the triangle
A = det([x1 y1 1
         x2 y2 1
         x3 y3 1])/2.0;

% find number of Gauss points
NP = length(xw(:,1)); 

z = 0.0;
for j = 1:NP
   L2 = xw(j,1);   
   L3 = xw(j,2);
   L1 = 1 - L2 - L3;
   
   x = x1*L1 + x2*L2 + x3*L3;
   y = y1*L1 + y2*L2 + y3*L3;
  
   z = z + f(x,y)*xw(j,3);
end

z = A*z;

return
