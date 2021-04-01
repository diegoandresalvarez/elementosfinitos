clear, clc

syms xi eta x1 x2 x3 x4 y3 y1 y2 y3 y4
N1 = (1 - xi  )*(1 - eta  )/4;
N2 = (1 + xi  )*(1 - eta  )/4;
N3 = (1 + xi  )*(1 + eta  )/4;
N4 = (1 - xi  )*(1 + eta  )/4;

N5 = (1 - xi^2)*(1 - eta  )/2;
N6 = (1 + xi  )*(1 - eta^2)/2;
N7 = (1 - xi^2)*(1 + eta  )/2;
N8 = (1 - xi  )*(1 - eta^2)/2;

%% Jacobian and inverse Jacobian
% isoparametric interpolation
x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

J = [ diff(x,xi)   diff(y,xi)
      diff(x,eta)  diff(y,eta) ];
disp('J = '); pretty(J)

invJ = inv(J);
% invJ = [ diff(xi,x)  diff(eta,x)
%          diff(xi,y)  diff(eta,y) ];

disp('inv(J) = '); disp(simplify(invJ))

%% Equation 12
syms x y
% Parameters of the shape funcions of a triangle
syms w1 w2 w3 w4 bx1 bx2 bx3 bx4 by1 by2 by3 by4

clear
%% Equation 27
syms xi eta phi5 phi6 phi7 phi8 dbs5 dbs6 dbs7 dbs8 L5 L6 L7 L8

dbn = [ dbs5 dbs6 dbs7 dbs8 ].';

gxiz5  = -L5*phi5*dbs5/6;
getaz6 = -L6*phi6*dbs6/6;
gxiz7  = +L7*phi7*dbs7/6;
getaz8 = +L8*phi8*dbs8/6;

gxiz  = (1-eta)*gxiz5  + (1+eta)*gxiz7;
getaz = (1-xi )*getaz8 + (1+xi )*getaz6;

syms j11 j12 j21 j22
gbar = [ j11 j12; j21 j22 ]*[ gxiz; getaz ];

bsdb = equationsToMatrix(gbar, dbn)

disp('bsdb = (1/6) * '); pretty(6*bsdb)
return

%%
clear
syms w1 w2 w3 w4 bx1 bx2 bx3 bx4 by1 by2 by3 by4 L5 L6 L7 L8 
syms phi4 phi5 phi6 phi7 phi8
syms x43 x32 x21 x14 y43 y32 y21 y14
c5 = x21/L5;      s5 = y21/L5;
c6 = x32/L6;      s6 = y32/L6;
c7 = x43/L7;      s7 = y43/L7;
c8 = x14/L8;      s8 = y14/L8;

Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 w4 bx4 by4 ].';

% eq 40b
% dbsk = (wj - wi + Lk*(ck*bxi + sk*byi)/2 + Lk*(ck*bxj + sk*byj)/2)/(-2*Lk*(1+phik)/3);
dbs5 = (w2 - w1 + L5*(c5*bx1 + s5*by1)/2 + L5*(c5*bx2 + s5*by2)/2)/(-2*L5*(1+phi5)/3);
dbs6 = (w3 - w2 + L6*(c6*bx2 + s6*by2)/2 + L6*(c6*bx3 + s6*by3)/2)/(-2*L6*(1+phi6)/3);
dbs7 = (w4 - w3 + L7*(c7*bx3 + s7*by3)/2 + L7*(c7*bx4 + s7*by4)/2)/(-2*L7*(1+phi7)/3);
dbs8 = (w1 - w4 + L8*(c8*bx4 + s8*by4)/2 + L8*(c8*bx1 + s8*by1)/2)/(-2*L8*(1+phi8)/3);

% Equation 41
An = equationsToMatrix([ dbs5; dbs6; dbs7; dbs8 ], Un);

% Equation 43
Adb = diag((2/3)*[ L5*(1+phi5) L6*(1+phi6) L7*(1+phi7) L8*(1+phi8) ])

% From equation 42, we obtain equation 44
Aw = Adb*An

%%
clear
syms gxz1 gyz1 gxz2 gyz2 gxz3 gyz3 gxz4 gyz4 gsz5 gsz6 gsz7 gsz8
syms C5 C6 C7 C8 S5 S6 S7 S8
syms A1 A2 A3 A4

syms xi eta x1 x2 x3 x4 y3 y1 y2 y3 y4
N1 = (1 - xi)*(1 - eta)/4;
N2 = (1 + xi)*(1 - eta)/4;
N3 = (1 + xi)*(1 + eta)/4;
N4 = (1 - xi)*(1 + eta)/4;

gxzbar = N1*gxz1 + N2*gxz2 + N3*gxz3 + N4*gxz4;
gyzbar = N1*gyz1 + N2*gyz2 + N3*gyz3 + N4*gyz4;
M6 = equationsToMatrix([ gxzbar; gyzbar ], [ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3; gxz4; gyz4 ]);

gxz1 = ( S8*gsz5 - S5*gsz8)/A1;
gyz1 = (-C8*gsz5 + S5*gsz8)/A1;

gxz2 = ( S5*gsz6 - S6*gsz5)/A2;
gyz2 = (-C5*gsz6 + S6*gsz5)/A2;

gxz3 = ( S6*gsz7 - S7*gsz6)/A3;
gyz3 = (-C6*gsz7 + S7*gsz6)/A3;

gxz4 = ( S7*gsz8 - S8*gsz7)/A4;
gyz4 = (-C7*gsz8 + S8*gsz7)/A4;
M5 = equationsToMatrix([ gxz1; gyz1; gxz2; gyz2; gxz3; gyz3; gxz4; gyz4 ], [ gsz5; gsz6; gsz7; gsz8 ])

syms phi5 phi6 phi7 phi8
M8 = diag(-(2/3)*[ phi5 phi6 phi7 phi8 ]);

Ngamma = M6*M5;
Bsdb = Ngamma*M8



%%
