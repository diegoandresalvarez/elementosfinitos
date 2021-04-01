clear, clc

syms xi eta x1 x2 x3 y1 y2 y3
lambda = 1 - xi - eta;
N1 = lambda;
N2 = xi;
N3 = eta;

%% Jacobian and inverse Jacobian
% isoparametric interpolation
x = N1*x1 + N2*x2 + N3*x3;
y = N1*y1 + N2*y2 + N3*y3;

J = [ diff(x,xi)   diff(y,xi)
      diff(x,eta)  diff(y,eta) ];
disp('J = '); pretty(J)

invJ = inv(J);
% invJ = [ diff(xi,x)  diff(eta,x)
%          diff(xi,y)  diff(eta,y) ];

% triangle area
Area = det([ 1 x1 y1          % Area del triangulo con vertices
             1 x2 y2          % (x1,y1), (x2,y2) y (x3,y3) numerados en el
             1 x3 y3])/2;     % sentido horario de las manecillas del reloj

disp('inv(J) = 1/(2*A) *'); pretty(simplify(invJ * 2*Area))

%% Equation 12
clear x y N1 N2 N3
% Parameters of the shape funcions of a triangle
syms w1 w2 w3 bx1 bx2 bx3 by1 by2 by3
syms x32 x13 x21 y32 y13 y21 A x y
a1 = x2*y3 - x3*y2;  b1 = -y32;  c1 = x32;
a2 = x3*y1 - x1*y3;  b2 = -y13;  c2 = x13;
a3 = x1*y2 - x2*y1;  b3 = -y21;  c3 = x21;

% Shape functions for a triangle of 3 nodes
N1 = (a1 + b1*x + c1*y)/(2*A);
N2 = (a2 + b2*x + c2*y)/(2*A);
N3 = (a3 + b3*x + c3*y)/(2*A);

%% Equation 12:
bx_13 = N1*bx1 + N2*bx2 + N3*bx3;
by_13 = N1*by1 + N2*by2 + N3*by3;

Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 ].';

Bbbeta = equationsToMatrix([ diff(bx_13,x)
                             diff(by_13,y)
                             diff(bx_13,y) + diff(by_13,x) ], Un);
disp('Bbbeta = 1/(2*A) * ');
pretty(Bbbeta * 2*A);

%% Equation 13
syms L4 L5 L6 cos4 cos5 cos6 sin4 sin5 sin6
x21 = cos4*L4;      y21 = sin4*L4;
x32 = cos5*L5;      y32 = sin5*L5;
x13 = cos6*L6;      y13 = sin6*L6;

P4 = 4*lambda*xi;
P5 = 4*xi*eta;
P6 = 4*lambda*eta;
%syms P4(xi,eta) P5(xi,eta) P6(xi,eta)
syms dbs4 dbs5 dbs6
bx_46 = P4*cos4*dbs4 + P5*cos5*dbs5 + P6*cos6*dbs6;
by_46 = P4*sin4*dbs4 + P5*sin5*dbs5 + P6*sin6*dbs6;

dbn = [ dbs4 dbs5 dbs6 ].';

dxi_dx = -y13/(2*A);    deta_dx = -y21/(2*A);
dxi_dy =  x13/(2*A);    deta_dy =  x21/(2*A);
diff_bx_p2_dx = diff(bx_46,xi)*dxi_dx + diff(bx_46,eta)*deta_dx;
diff_bx_p2_dy = diff(bx_46,xi)*dxi_dy + diff(bx_46,eta)*deta_dy;
diff_by_p2_dx = diff(by_46,xi)*dxi_dx + diff(by_46,eta)*deta_dx;
diff_by_p2_dy = diff(by_46,xi)*dxi_dy + diff(by_46,eta)*deta_dy;

Bbdbeta = equationsToMatrix([ diff_bx_p2_dx
                              diff_by_p2_dy
                              diff_bx_p2_dy + diff_by_p2_dx ], dbn);
                         
disp('Bbdbeta = 1/(2*A) * ');
disp(simplify(expand(Bbdbeta*2*A)))

%%
clear
syms w1 w2 w3 bx1 bx2 bx3 by1 by2 by3 L4 L5 L6 c4 c5 c6 s4 s5 s6 phi4 phi5 phi6
syms x32 x13 x21 y32 y13 y21 
c4 = x21/L4;      s4 = y21/L4;
c5 = x32/L5;      s5 = y32/L5;
c6 = x13/L6;      s6 = y13/L6;

Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 ].';

% eq 40b
% dbsk = (wj - wi + Lk*(ck*bxi + sk*byi)/2 + Lk*(ck*bxj + sk*byj)/2)/(-2*Lk*(1+phik)/3);
dbs4 = (w2 - w1 + L4*(c4*bx1 + s4*by1)/2 + L4*(c4*bx2 + s4*by2)/2)/(-2*L4*(1+phi4)/3);
dbs5 = (w3 - w2 + L5*(c5*bx2 + s5*by2)/2 + L5*(c5*bx3 + s5*by3)/2)/(-2*L5*(1+phi5)/3);
dbs6 = (w1 - w3 + L6*(c6*bx3 + s6*by3)/2 + L6*(c6*bx1 + s6*by1)/2)/(-2*L6*(1+phi6)/3);

% Equation 41
An = equationsToMatrix([ dbs4; dbs5; dbs6 ], Un);

% Equation 43
Adb = diag((2/3)*[ L4*(1+phi4) L5*(1+phi5) L6*(1+phi6) ])

% From equation 42, we obtain equation 44
Aw = Adb*An

%%
