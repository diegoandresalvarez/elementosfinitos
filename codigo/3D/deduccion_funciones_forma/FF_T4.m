clear, clc
%% Programa para verificar las funciones de forma del elemento tetraedrico 
%% de cuatro nodos que aparecen en la literatura

%% Se definen las variables simbolicas
syms x y z Vol
syms x1 x2 x3 x4 
syms y1 y2 y3 y4 
syms z1 z2 z3 z4 
syms u1 u2 u3 u4 
syms alpha_1 alpha_2 alpha_3 alpha_4

%% Resuelve el sistema de ecuaciones por alpha_1, alpha_2, alpha_3 y alpha_4
r = solve(u1 == alpha_1 + alpha_2*x1 + alpha_3*y1 + alpha_4*z1, ...
          u2 == alpha_1 + alpha_2*x2 + alpha_3*y2 + alpha_4*z2, ...
          u3 == alpha_1 + alpha_2*x3 + alpha_3*y3 + alpha_4*z3, ...
          u4 == alpha_1 + alpha_2*x4 + alpha_3*y4 + alpha_4*z4, ...
          alpha_1, alpha_2, alpha_3, alpha_4);

%% Establece de nuevo u
u = r.alpha_1 + r.alpha_2*x + r.alpha_3*y + r.alpha_4*z;

%% Se separa u en su numerador y su denominador
[num,den] = numden(u);
num = -num;
den = -den;

%% Nosotros sabemos que el denominador es seis veces el volumen del tetraedro
% A continuacion se verificara:
V = det([ 1 x1 y1 z1          % Volumen del tetraedro con vertices
          1 x2 y2 z2          % (x1,y1,z1), ..., (x4,y4,z4) 
          1 x3 y3 z3          % suponiendo que (x1,y1,z1), ..., (x3,y3,z3)          
          1 x4 y4 z4 ])/6;    % se numeraron en sentido antihorario cuando se mira desde (x4,y4,z4)
VV = matlabFunction(V, 'Vars', [x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4]);
      
disp('El 0 en el residuo significa que den == 6*V');
disp(simplify(den - 6*V))

%% Se vuelve a escribir u, pero con el denominador expresado como 6*Vol
% ya que en el paso anterior establecimos que den es igual a 6*V
u = num/(6*Vol);

%% Se factoriza u1, u2, u3 y u4
u = collect(u,[u1, u2, u3, u4]);

%% Se verifican finalmente las formulas
a1 =   det([x2 y2 z2; x3 y3 z3; x4 y4 z4]);
a2 =  -det([x3 y3 z3; x4 y4 z4; x1 y1 z1]);
a3 =   det([x4 y4 z4; x1 y1 z1; x2 y2 z2]);
a4 =  -det([x1 y1 z1; x2 y2 z2; x3 y3 z3]);

b1 = -det([ 1 y2 z2;  1 y3 z3;  1 y4 z4]); % OJO CON LOS SIGNOS
b2 =  det([ 1 y3 z3;  1 y4 z4;  1 y1 z1]); % NO DAN LO MISMO QUE LA 
b3 = -det([ 1 y4 z4;  1 y1 z1;  1 y2 z2]); % LITERATURA
b4 =  det([ 1 y1 z1;  1 y2 z2;  1 y3 z3]);

c1 = -det([x2  1 z2; x3  1 z3; x4  1 z4]);
c2 =  det([x3  1 z3; x4  1 z4; x1  1 z1]);
c3 = -det([x4  1 z4; x1  1 z1; x2  1 z2]);
c4 =  det([x1  1 z1; x2  1 z2; x3  1 z3]);

d1 = -det([x2 y2  1; x3 y3  1; x4 y4  1]);
d2 =  det([x3 y3  1; x4 y4  1; x1 y1  1]);
d3 = -det([x4 y4  1; x1 y1  1; x2 y2  1]);
d4 =  det([x1 y1  1; x2 y2  1; x3 y3  1]);

%{
% Formulas incorrectas que aparecen en los libros
a1 =   det([x2 y2 z2; x3 y3 z3; x4 y4 z4]);
a2 =   det([x3 y3 z3; x4 y4 z4; x1 y1 z1]);
a3 =   det([x4 y4 z4; x1 y1 z1; x2 y2 z2]);
a4 =   det([x1 y1 z1; x2 y2 z2; x3 y3 z3]);

b1 = -det([ 1 y2 z2;  1 y3 z3;  1 y4 z4]); % OJO CON LOS SIGNOS
b2 = -det([ 1 y3 z3;  1 y4 z4;  1 y1 z1]); % NO DAN LO MISMO QUE LA 
b3 = -det([ 1 y4 z4;  1 y1 z1;  1 y2 z2]); % LITERATURA
b4 = -det([ 1 y1 z1;  1 y2 z2;  1 y3 z3]);

c1 = -det([x2  1 z2; x3  1 z3; x4  1 z4]);
c2 = -det([x3  1 z3; x4  1 z4; x1  1 z1]);
c3 = -det([x4  1 z4; x1  1 z1; x2  1 z2]);
c4 = -det([x1  1 z1; x2  1 z2; x3  1 z3]);

d1 = -det([x2 y2  1; x3 y3  1; x4 y4  1]);
d2 = -det([x3 y3  1; x4 y4  1; x1 y1  1]);
d3 = -det([x4 y4  1; x1 y1  1; x2 y2  1]);
d4 = -det([x1 y1  1; x2 y2  1; x3 y3  1]);
%}

%% Se arman las funciones de forma
N  = cell(4,1);
N{1} = (a1 + b1*x + c1*y + d1*z)/(6*Vol);
N{2} = (a2 + b2*x + c2*y + d2*z)/(6*Vol);
N{3} = (a3 + b3*x + c3*y + d3*z)/(6*Vol);
N{4} = (a4 + b4*x + c4*y + d4*z)/(6*Vol);

%% Se arma la funcion de desplazamiento
uu = N{1}*u1 + N{2}*u2 + N{3}*u3 + N{4}*u4;

disp('El 0 en el residuo significa que ambas formulaciones coinciden');
disp(simplify(u-uu));

%% Evaluacion de las funciones de forma en los vertices
xnod = [ ...
   -1.0000   -1.0000         0
   -0.9000   -0.9000         0
   -0.8000   -1.0000         0
   -1.0000   -0.8000         0
   -1.0000   -0.9000    0.1000
   -0.9000   -1.0000    0.1000
   -1.0000   -1.0000    0.2000
   -0.8000   -1.0000    0.2000
   -0.9000   -0.9000    0.2000 ];

i = 1; x1 = xnod(i,1); y1 = xnod(i,2); z1 = xnod(i,3);
i = 5; x2 = xnod(i,1); y2 = xnod(i,2); z2 = xnod(i,3);
i = 7; x3 = xnod(i,1); y3 = xnod(i,2); z3 = xnod(i,3);
i = 6; x4 = xnod(i,1); y4 = xnod(i,2); z4 = xnod(i,3);

% Vertices para ensayar (se esta ensayando actualmente la cuarta fila)
%     1     4     5     2
%     1     5     6     2
%     1     6     3     2
%     1     5     7     6

%% Se calcula el volumen del tetrahedro
fprintf('Vol = %d\n\n',VV(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4));

%% Se evaluan todas las funciones de forma en cada uno de sus vertices
N  = cell(4,1);
N{1} = (a1 + b1*x + c1*y + d1*z)/(6*V);
N{2} = (a2 + b2*x + c2*y + d2*z)/(6*V);
N{3} = (a3 + b3*x + c3*y + d3*z)/(6*V);
N{4} = (a4 + b4*x + c4*y + d4*z)/(6*V);

for i = 1:4
   NN = matlabFunction(N{i}, 'Vars', {'x','y','z', 'x1','y1','z1', 'x2','y2','z2', 'x3','y3','z3', 'x4','y4','z4'});
    
   [ NN(x1,y1,z1, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)
     NN(x2,y2,z2, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)
     NN(x3,y3,z3, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)
     NN(x4,y4,z4, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4) ]
end

%% Se verifica la condición de cuerpo rígido: sum(N) == 1
suma = 0;
for i = 1:4
   suma = suma + N{i};
end
fprintf('\nSe verifica la condición de cuerpo rígido: sum(N) == ');
disp(simplify(suma));
