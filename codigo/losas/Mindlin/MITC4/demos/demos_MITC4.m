%% Demonstration of the equations relative to the MITC4 finite element in:
%
% Katili, I., Bato, J.-L., Maknun, J. and Lardeur, P. (2018), A comparative 
% formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements with new 
% numerical results based on s-norm tests. Computers & Structures, 204:
% 48-64. https://doi.org/10.1016/j.compstruc.2018.04.001

%%
clear, clc
disp('EQUATIONS Bs KATILI ET. AL.')
disp('*** Calculation of the inverse Jacobian j = inv(J) ***')

% Jacobian and inverse Jacobian
% isoparametric interpolation

syms xi eta x1 x2 x3 x4 y1 y2 y3 y4
N1 = (1 - xi  )*(1 - eta  )/4; % = h3 Bathe
N2 = (1 + xi  )*(1 - eta  )/4; % = h4 Bathe
N3 = (1 + xi  )*(1 + eta  )/4; % = h1 Bathe
N4 = (1 - xi  )*(1 + eta  )/4; % = h2 Bathe

x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

disp('J = '); 
J = [ diff(x,xi)   diff(y,xi)
      diff(x,eta)  diff(y,eta) ]

disp('det(J) = ');
detJ = det(J)

invJ = inv(J);
% invJ = [ diff(xi,x)  diff(eta,x)
%          diff(xi,y)  diff(eta,y) ];

disp('4*inv(J)*det(J) = '); collect(4*invJ*detJ,{'r','s','xi','eta'})

%% Equation Katili, equation 44
disp('*** DEMO KATILI, EQUATION 44 ***')

syms gxiz5 getaz6 gxiz7 getaz8
g_xieta = [ gxiz5 getaz6 gxiz7 getaz8 ].';

% Bathe, equation 6
% grz = (1/2)*(1-s)*grzC + (1/2)*(1+s)*grzA;
% gsz = (1/2)*(1-r)*gszB + (1/2)*(1+r)*gszD;

gxiz  = (1/2)*(1-eta)*gxiz5  + (1/2)*(1+eta)*gxiz7;
getaz = (1/2)*(1-xi )*getaz8 + (1/2)*(1+xi )*getaz6;

gbar = [ gxiz; getaz ];

Ng = equationsToMatrix(gbar, g_xieta);

disp('Ngamma = (1/2)*'); disp(2*Ng)

%%
disp('*** DEMO KATILI, EQUATION 50 ***')

syms w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 w4 bx4 by4
syms x43 x32 x21 x14 y43 y32 y21 y14 L5 L6 L7 L8 
c5 = x21/L5;      s5 = y21/L5;
c6 = x32/L6;      s6 = y32/L6;
c7 = x43/L7;      s7 = y43/L7;
c8 = x14/L8;      s8 = y14/L8;

Un = [ w1 bx1 by1 w2 bx2 by2 w3 bx3 by3 w4 bx4 by4 ].';

% eq 49
% gsk = (wj - wi)/Lk + (ck*bxi + sk*byi)/2 + (ck*bxj + sk*byj)/2;
gs5   = (w2 - w1)/L5 + (c5*bx1 + s5*by1)/2 + (c5*bx2 + s5*by2)/2;
gs6   = (w3 - w2)/L6 + (c6*bx2 + s6*by2)/2 + (c6*bx3 + s6*by3)/2;
gs7   = (w4 - w3)/L7 + (c7*bx3 + s7*by3)/2 + (c7*bx4 + s7*by4)/2;
gs8   = (w1 - w4)/L8 + (c8*bx4 + s8*by4)/2 + (c8*bx1 + s8*by1)/2;

% Equation 34/50
Au = equationsToMatrix([ gs5; gs6; gs7; gs8 ], Un);
disp('Au = (1/2)*'); pretty(2*Au);

%%
disp('*** DEMO KATILI, EQUATION 52 ***')

% Equation 46
Ag = diag([ L5 L6 -L7 -L8 ]/2)

% Equation 52
Ng_Ag_Au = simplify(expand(Ng*Ag*Au));
disp('(Ng*Ag*Au)^T = '); Ng_Ag_Au.'

Bs = invJ*Ng_Ag_Au


%%
disp('************************************************************************')
disp('EQUATIONS Bs BATHE-DVORKIN')
clear
syms r s
syms w1 tx1 ty1 w2 tx2 ty2 w3 tx3 ty3 w4 tx4 ty4 
syms x14 x21 x32 x43 y14 y21 y32 y43 

%{
eq_grz = (1+s)*((w1 - w2)/2 + ((x1 - x2)/4)*(ty1 + ty2) - ((y1 - y2)/4)*(tx1 + tx2)) + ...
         (1-s)*((w4 - w3)/2 + ((x4 - x3)/4)*(ty4 + ty3) - ((y4 - y3)/4)*(tx4 + tx3));
eq_gsz = (1+r)*((w1 - w4)/2 + ((x1 - x4)/4)*(ty1 + ty4) - ((y1 - y4)/4)*(tx1 + tx4)) + ...
         (1-r)*((w2 - w3)/2 + ((x2 - x3)/4)*(ty2 + ty3) - ((y2 - y3)/4)*(tx2 + tx3));
%}
eq_grz = (1+s)*((w1 - w2)/2 - (x21/4)*(ty1 + ty2) + (y21/4)*(tx1 + tx2)) + ...
         (1-s)*((w4 - w3)/2 + (x43/4)*(ty4 + ty3) - (y43/4)*(tx4 + tx3));
eq_gsz = (1+r)*((w1 - w4)/2 + (x14/4)*(ty1 + ty4) - (y14/4)*(tx1 + tx4)) + ...
         (1-r)*((w2 - w3)/2 - (x32/4)*(ty2 + ty3) + (y32/4)*(tx2 + tx3));

ae = [ w1 tx1 ty1 w2 tx2 ty2 w3 tx3 ty3 w4 tx4 ty4 ].';
Bgrz = simplify(expand(equationsToMatrix(eq_grz, ae))).'
Bgsz = simplify(expand(equationsToMatrix(eq_gsz, ae))).'

