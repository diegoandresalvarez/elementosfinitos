% Ejemplo 2.4. Onate, tomo II

clear  % borrar memoria y pantalla

%% defino las constantes y variables
Y  = 1; TH = 2; 
T1 = 1; T2 = 2; EB = 3;

syms b h E I G Aast L P lambda nu  % define las variables simbolicas

% Cambie el numero de EFs para generar las ecuaciones
nef = 8;                           % numero de elementos finitos (EF)
nno = nef + 1;                     % numero de nodos
ngdl = 2*nno;                      % numero de grados de libertad

% OJO: Onate programo lo siguiente:
Le = L;                            % longitud de la barra
% Yo hubiera colocado Le = L/nef. Cuando lo hago asi, el codigo no funciona
% para compensar, la altura se programa como nef*h

gdl = [ (1:2:ngdl)'  (2:2:ngdl)']; % grados de libertad
LaG = [ (1:(nno-1))' (2:nno)'   ]; % definicion de EFs con respecto a nodos

%% Relaciono apoyos
%  gdl           desplazamiento
apoyos = [ ...
   gdl(  1,Y)    0    % m
   gdl(  1,TH)   0 ]; % rad   

%% grados de libertad del desplazamiento conocidos y desconocidos
c  = apoyos(:, 1);           % GDL conocidos
d =  setdiff((1:ngdl)', c);  % GDL desconocidos
ac = sym(apoyos(:, 2));      % desplazamientos conocidos

%% Relacion de cargas puntuales
f = sym(zeros(ngdl,1));   % vector de fuerzas nodales equivalentes global
f(gdl(nno, Y)) = P;       % carga puntual en x = L

%% se separa la memoria para las matrices de rigidez
KT1 = sym(zeros(ngdl));
KT2 = sym(zeros(ngdl));
KEB = sym(zeros(ngdl));

%% ensamblo la matriz de rigidez global (K mayuscula)
for e = 1:nef         % para cada una de las barras e = 1, 2 y 3
   idx = [ gdl(LaG(e,1), :) gdl(LaG(e,2), :) ];
   
   Kb_T = (E*I/L) * [ 0,  0, 0,  0
                      0,  1, 0, -1
                      0,  0, 0,  0
                      0, -1, 0,  1 ];
   
   Ks1_T = (G*Aast/Le) * [   1,   Le/2,    -1,   Le/2
                          Le/2, Le^2/4, -Le/2, Le^2/4
                            -1,  -Le/2,     1,  -Le/2
                          Le/2, Le^2/4,  -Le/2, Le^2/4 ];

   Ks2_T = (G*Aast/Le) * [   1,   Le/2,    -1,    Le/2
                          Le/2, Le^2/3, -Le/2,  Le^2/6
                            -1,  -Le/2,     1,   -Le/2
                          Le/2, Le^2/6,  -Le/2, Le^2/3 ];

   K_EB = (E*I/L^3) * [  12,   6*L,  -12,   6*L
                        6*L, 4*L^2, -6*L, 2*L^2
                        -12,  -6*L,   12,  -6*L
                        6*L, 2*L^2, -6*L, 4*L^2 ];

   KT1(idx,idx) = KT1(idx,idx) + Kb_T + Ks1_T;
   KT2(idx,idx) = KT2(idx,idx) + Kb_T + Ks2_T;
   KEB(idx,idx) = KEB(idx,idx) + K_EB;                          
end

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |     recuerde siempre que qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |     en este caso en particular fd=0

a = cell(1,3);
q = cell(1,3);
for i = 1:3
   switch i
      case 1
         K = KT1;
      case 2
         K = KT2;
      case 3
         K = KEB;         
   end

   Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
   Kdc = K(d,c); Kdd = K(d,d); fc = f(d);   

   %% resuelvo el sistema de ecuaciones
   % recuerde que \ es para resolver el sistema de ecuaciones eficientemente
   ad = simplify(Kdd\(fc - Kdc*ac));  % = linsolve(Kdd, fc-Kdc*ac)
   qd = simplify(Kcc*ac + Kcd*ad);

   %% formo los vectores de desplazamientos (a) y fuerzas (q)
   aa = sym(zeros(ngdl,1)); qq = sym(zeros(ngdl,1));  % separo la memoria
   aa(c) = ac;              qq(c) = qd;
   aa(d) = ad;            % qq(d) = qc = 0; 
   
   a{i} = aa;
   q{i} = qq;
end

%%
h = nef*L/lambda;

rwT1 = a{T1}(end-1)/a{EB}(end-1);
rwT1 = simplify(subs(rwT1, {I, Aast, G}, {b*h^3/12, 5*b*h/6, E/(2*(1+nu))}))
limit(rwT1, lambda, inf)
%%
rwT2 = a{T2}(end-1)/a{EB}(end-1);
rwT2 = simplify(subs(rwT2, {I, Aast, G}, {b*h^3/12, 5*b*h/6, E/(2*(1+nu))}))
limit(rwT2, lambda, inf)

return

%% A partir de los resultados se hizo el siguiente programa:
clear, clc, close all

lambda = linspace(0, 20, 5000);
nu = 0.25;

rwT1_1 = (15*lambda.^2 + 12*nu + 12)./(20*lambda.^2);
rwT2_1 = (12*(nu + 1)*(5*lambda.^2 + 3*nu + 3))./(5*lambda.^2.*(5*lambda.^2 + 12*nu + 12));
rwT1_2 = (75*lambda.^2 + 48*nu + 48)./(80*lambda.^2);
rwT2_2 = (48*(nu + 1)*(5*lambda.^2 + 3*nu + 3))./(5*lambda.^2.*(5*lambda.^2 + 48*nu + 48));
rwT1_4 = (315*lambda.^2 + 192*nu + 192)./(320*lambda.^2);
rwT2_4 = (192*(nu + 1)*(5*lambda.^2 + 3*nu + 3))./(5*lambda.^2.*(5*lambda.^2 + 192*nu + 192));
rwT1_8 = (1275*lambda.^2 + 768*nu + 768)./(1280*lambda.^2);
rwT2_8 = (768*(nu + 1)*(5*lambda.^2 + 3*nu + 3))./(5*lambda.^2.*(5*lambda.^2 + 768*nu + 768));

figure
hold on;
plot(lambda, rwT1_1, 'r', lambda, rwT2_1, 'r--');
plot(lambda, rwT1_2, 'c', lambda, rwT2_2, 'c--');
plot(lambda, rwT1_4, 'm', lambda, rwT2_4, 'm--');
plot(lambda, rwT1_8, 'b', lambda, rwT2_8, 'b--');
legend('Ks (1 point), 1 EF', 'Ks (2 point), 1 EF', ...
       'Ks (1 point), 2 EF', 'Ks (2 point), 2 EF', ...
       'Ks (1 point), 4 EF', 'Ks (2 point), 4 EF', ...  
       'Ks (1 point), 8 EF', 'Ks (2 point), 8 EF');
ylim([0, 4]);
grid on;

xlabel('$$\lambda$$','interpreter','latex')
ylabel('$$r_w$$','interpreter','latex')

%% Y calculamos los limites
clear
nu = 0.25;
syms lambda

rwT1_1 = (15*lambda^2 + 12*nu + 12)/(20*lambda^2);
rwT2_1 = (12*(nu + 1)*(5*lambda^2 + 3*nu + 3))/(5*lambda^2*(5*lambda^2 + 12*nu + 12));
rwT1_2 = (75*lambda^2 + 48*nu + 48)/(80*lambda^2);
rwT2_2 = (48*(nu + 1)*(5*lambda^2 + 3*nu + 3))/(5*lambda^2*(5*lambda^2 + 48*nu + 48));
rwT1_4 = (315*lambda^2 + 192*nu + 192)/(320*lambda^2);
rwT2_4 = (192*(nu + 1)*(5*lambda^2 + 3*nu + 3))/(5*lambda^2*(5*lambda^2 + 192*nu + 192));
rwT1_8 = (1275*lambda^2 + 768*nu + 768)/(1280*lambda^2);
rwT2_8 = (768*(nu + 1)*(5*lambda^2 + 3*nu + 3))/(5*lambda^2*(5*lambda^2 + 768*nu + 768));

limit(rwT1_1, lambda, inf)   % = 3/4     = 0.7500
limit(rwT2_1, lambda, inf)   % = 0
limit(rwT1_2, lambda, inf)   % = 15/16   = 0.9375
limit(rwT2_2, lambda, inf)   % = 0
limit(rwT1_4, lambda, inf)   % = 63/64   = 0.9844
limit(rwT2_4, lambda, inf)   % = 0
limit(rwT1_8, lambda, inf)   % = 255/256 = 0.9961
limit(rwT2_8, lambda, inf)   % = 0
