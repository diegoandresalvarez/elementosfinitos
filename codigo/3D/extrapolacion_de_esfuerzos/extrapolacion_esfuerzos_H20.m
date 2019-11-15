clear, clc

%% Numeracion local del EF serendipito hexaedrico de 20 nodos:
%
%      ^ eta        para zeta = +1
%      |
%      |
% 19--18--17                                             ^ zeta
%  |   |   |                                             |
% 20---+--16----> xi                                     |            
%  |   |   |                                       13---14----15   zeta = +1
% 13--14--15                                      /|    /     /|
%                                                20----+----16 | 
%                                               /  |  /     /  | 
%      ^ eta        para zeta = 0              19---18----17   |      
%      |                                       |   |       |   |
%      |                                       |   09----+-|--10   zeta =  0
% 12---+--11                                   |  /|    /  |  /|
%  |   |   |                                   | +-----+---|-+ | 
%  +---+---+----> xi                           |/  |  /    |/  |
%  |   |   |                                   12----+----11   |      
%  9---+--10                                   |   |       |   | 
%                                              |   01---02-|--03   zeta = -1
%                                              |  /     /  |  /
% Numeracion local:                            | 08----+---|04-------> xi
%      ^ eta        para zeta = -1             |/     /    |/
%      |                                       07---06----05         
%      |                                            /
%  7---6---5                                       /
%  |   |   |                                      / eta
%  8---+---4----> xi
%  |   |   |
%  1---2---3

%% Coordenadas de los nodos:
nod = [ ...
%  xi   eta  zeta   nodo
   -1   -1   -1   %  1
    0   -1   -1   %  2
    1   -1   -1   %  3
    1    0   -1   %  4
    1    1   -1   %  5
    0    1   -1   %  6
   -1    1   -1   %  7
   -1    0   -1   %  8
   -1   -1    0   %  9
    1   -1    0   % 10
    1    1    0   % 11
   -1    1    0   % 12
   -1   -1    1   % 13
    0   -1    1   % 14
    1   -1    1   % 15
    1    0    1   % 16
    1    1    1   % 17
    0    1    1   % 18
   -1    1    1   % 19
   -1    0    1 ];% 20

%% Se define la cuadratura de Gauss Legendre a utilizar
[x_gl, w_gl] = gausslegendre_quad_hexa(2);
x_gl = sym(x_gl);
n_gl = length(w_gl);

%% Se crea la matriz A1
xi   = x_gl(:,1); 
eta  = x_gl(:,2);
zeta = x_gl(:,3);
A1 = [ ones( 8,1) xi eta zeta xi.*eta xi.*zeta eta.*zeta xi.*eta.*zeta ];

%% Se crea la matriz A2
xi   = nod(:,1);
eta  = nod(:,2);
zeta = nod(:,3);
A2 = [ ones(20,1) xi eta zeta xi.*eta xi.*zeta eta.*zeta xi.*eta.*zeta ];

%% Se reporta la matriz de extrapolacion
fprintf('La matriz de extrapolacion del EF H20 es:\n')
A = simplify(A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))