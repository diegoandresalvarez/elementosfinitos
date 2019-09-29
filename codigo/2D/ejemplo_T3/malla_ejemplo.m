xnod = [ 0     0    % posicion de los nodos
         0.2   0
         0.4   0    % la fila representa el numero del nodo
         0.6   0    % la columna representa la coordenada X=1 o Y=2
         0.8   0
         0     0.14
         0.2   0.14 
         0.4   0.14
         0.6   0.14
         0.8   0.14
         0     0.28
         0.2   0.28 
         0.4   0.28
         0.6   0.28
         0.8   0.28 ];

nno = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 2*nno;       % numero de grados de libertad (dos por nodo)
gdl = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad
 
LaG = [  1  2  7    % definicion de elementos finitos con respecto a nodos
         7  6  1
         6  7 12    % la fila representa el numero del elemento
        12 11  6    % la columna representa el numero de nodo local 1,2 o 3
         2  3  8
         8  7  2
         7  8 13
        13 12  7
         3  4  9
         9  8  3
         8  9 14
        14 13  8
         4  5 10
        10  9  4
         9 10 15
        15 14  9 ];
      
nef = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se definen las restricciones 
%             gdl       desplazamiento(m)
restric = [   gdl(1,X)          0
              gdl(1,Y)          0
              gdl(5,Y)          0  ];
           
%% Se definen las cargas distribuidas 
%             [ elemento lado tix  tiy  tjx  tjy  ]
carga_distr = [     4     12    0 -7.5    0 -10.0
                    8     12    0 -5.0    0  -7.5 ];
                 
nlcd = size(carga_distr,1); % numero de lados con carga distribuida
         
%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
f(gdl(14,X)) = 7000*cosd(60);  % carga puntual en el nodo 14 dir X
f(gdl(14,Y)) = 7000*sind(60);  % carga puntual en el nodo 14 dir Y


%{
% ESTE CODIGO SIRVE PARA CALCULAR "A PEDAL" LAS FUERZAS NODALES 
% SUPERFICIALES EQUIVALENTES

% (toca hacer "a pedal" este pedazo porque programarlo de modo mas
% automatico tomaria muchas lineas)
ft = sparse(ngdl,1); % fuerzas nodales equivalentes de cargas superficiales

syms x
carga = poly2sym(polyfit([0 0.4],[-10 -5],1), x); 
% esto es: carga = @(x) (25/2)kN/m^2 x - 10 kN/m (sobre borde 11-12-13)

% Recuerde que la integral de linea se debe hacer sobre el contorno del
% elemento. En este caso la funcion de forma bidimensional sobre el
% contorno se vuelve una funcion de forma similar a la funcion de forma 
% unidimensional que pasa por dos nodos similar a polyfit([0 0.2],[1 0],1)

% Elemento 4 11----12----13
ft(gdl(11,Y)) = ... % Nodo 2 local = 11 global
   double(int(poly2sym(polyfit([0 0.2],[1 0],1), x)*carga*te,x,0,0.2));

ft(gdl(12,Y)) = ... % Nodo 1 local = 12 global
   double(int(poly2sym(polyfit([0 0.2],[0 1],1), x)*carga*te,x,0,0.2));

% Elemento 8
ft(gdl(12,Y)) = ft(gdl(12,Y)) + ...   % Nodo 2 local = 12 global
    double(int(poly2sym(polyfit([0.2 0.4],[1 0],1), x)*carga*te,x,0.2,0.4));

ft(gdl(13,Y)) = ...                   % Nodo 1 local = 13 global
   double(int(poly2sym(polyfit([0.2 0.4],[0 1],1), x)*carga*te,x,0.2,0.4));
%}