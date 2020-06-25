E  = 6.825e7;        % modulo de elasticidad del solido (Pa)
nu = 0.3;            % coeficiente de Poisson
qdistr = -90;        % carga (N/m^2)

%% Definimos la geometria de la losa 
% Se genero con generacion_malla_cascaron_semiesfera_S8.m
r998  = load('semiesfera_curvo_S8_9_98', 'LaG','xnod');
r1002 = load('semiesfera_curvo_S8_10_02','LaG','xnod');
if ~isequal(r998.LaG, r1002.LaG)
   error('Las matrices LaG no coinciden');
end
LaG = r998.LaG;
xnod_lo = r998.xnod;
xnod_up = r1002.xnod;
xnod = (xnod_lo + xnod_up)/2;

% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nnoef = size(LaG,2);           % numero de nodos del EF

%% Se crea el vector tt
tt = zeros(5*nnoef,1);
tt(ww:5:end) = qdistr;

%% Se definen las restricciones 
% determino los grados de libertad correspondientes a los bordes
base = find(abs(xnod(:,Z) - 0) < 1e-4);

%           nodo direccion desplazamiento(m)
restricciones = [...
   base repmat([uu 0],size(base))
   base repmat([vv 0],size(base))
   base repmat([ww 0],size(base))
%   base repmat([t1 0],size(base))
%   base repmat([t2 0],size(base))
];

%% BYE BYE!!