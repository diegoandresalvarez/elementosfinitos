% Calculo de las funciones de forma del elemento placa MZC

clear, clc, close all

syms x y xi eta a b x0 y0 E nu t q

PT = [ 1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3 ];
PT_x = diff(PT,x);
PT_y = diff(PT,y);

A = [ ...
      subs(PT,   {x, y}, {-1, -1})
      subs(PT_x, {x, y}, {-1, -1}) % nodo 1
      subs(PT_y, {x, y}, {-1, -1})
      
      subs(PT,   {x, y}, { 1, -1})
      subs(PT_x, {x, y}, { 1, -1}) % nodo 2
      subs(PT_y, {x, y}, { 1, -1})
      
      subs(PT,   {x, y}, { 1,  1})
      subs(PT_x, {x, y}, { 1,  1}) % nodo 3
      subs(PT_y, {x, y}, { 1,  1})
      
      subs(PT,   {x, y}, {-1,  1})
      subs(PT_x, {x, y}, {-1,  1}) % nodo 4
      subs(PT_y, {x, y}, {-1,  1}) ];

N = PT/A;  % N = PT*inv(A);

% reescribo las funciones de forma en terminos de xi y eta
N = subs(N,{x, y},{xi, eta}); 

% hallo las derivadas
d3N_dxi3 = diff(N,xi,3)/(a^3); d3N_deta3 = diff(N,eta,3)/(b^3);
d3N_dxideta2 = diff(N,xi,eta,eta)/(a*b^2);
d3N_dxi2deta = diff(N,xi,xi,eta)/(a^2*b);

% reemplazo en la matriz R
R = [d3N_dxi3 + d3N_dxideta2
     d3N_dxi2deta + d3N_deta3];
R = simplify(R);

% Q = -D*R*a
R

