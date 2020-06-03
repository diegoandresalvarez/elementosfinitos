function r = lhs(expr)
% LHS Returns the left hand side an expression
% LHS(sym('[1-x^2==2*y'; 1+x^2==x+y]')) = [1-x^2; 1+x^2]

LHS = @(e)['lhs(',char(e),')'];
[m,n] = size(expr);
r = sym(zeros(m,n));
for i=1:m
    for j=1:n
        r(i,j) = evalin(symengine, LHS(expr(i,j)));
    end
end  
