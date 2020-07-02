E  = 6.825e7;        % modulo de elasticidad del solido (Pa)
nu = 0.3;            % coeficiente de Poisson
t  = 0.04;           % espesor del cascaron (m)
qdistr = -90;        % carga (N/m^2)

%% Definimos la geometria de la losa 
% Se genero con cascaron_scordelli_lo.m
load semiesfera_Q4 LaG xnod
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nnoef = size(LaG,2);           % numero de nodos del EF

%% Se crea el vector tt
tt = zeros(6*nnoef,1);
tt(ww:6:end) = qdistr;

%% Se definen las restricciones 
% determino los grados de libertad correspondientes a los bordes
base = find(abs(xnod(:,Z) - 0) < 1e-4);

%           nodo direccion desplazamiento(m)
restricciones = [...
   base repmat([uu 0],size(base))
   base repmat([vv 0],size(base))
   base repmat([ww 0],size(base))
   base repmat([tx 0],size(base))
   base repmat([ty 0],size(base))
   base repmat([tz 0],size(base))   
];

%% BYE BYE!!