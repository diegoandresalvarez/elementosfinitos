E  = 6.825e7;        % modulo de elasticidad del solido (Pa)
nu = 0.3;            % coeficiente de Poisson
t  = 0.3;           % espesor del cascaron (m)
qdistr = -90;        % carga (N/m^2)

%% Definimos la geometria de la losa 
% Se genero con cascaron_scordelli_lo.m
%load scordelli_lo_Q4_malla_15_10 LaG xnod  % converge
load scordelli_lo_Q4_malla_30_10 LaG xnod  % converge
%load scordelli_lo_Q4_malla_30_20 LaG xnod  % no converge (problemas de coplaneraidad -- lamina rebajada)


% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos

nnoef = size(LaG,2);           % numero de nodos del EF

%% Se crea el vector tt
tt = zeros(6*nnoef,1);
tt(ww:6:end) = qdistr;

%% Se definen las restricciones 
% determino los grados de libertad correspondientes a los bordes
lado_ypos = find(abs(xnod(:,Y) - (+25*sind(40))) < 1e-4);
lado_yneg = find(abs(xnod(:,Y) - (-25*sind(40))) < 1e-4);

%           nodo direccion desplazamiento(m)
restricciones = [...
   lado_ypos repmat([uu 0],size(lado_ypos))
   lado_ypos repmat([vv 0],size(lado_ypos))
   lado_ypos repmat([ww 0],size(lado_ypos))
   lado_yneg repmat([uu 0],size(lado_ypos))
   lado_yneg repmat([vv 0],size(lado_ypos))
   lado_yneg repmat([ww 0],size(lado_ypos))
];

%% BYE BYE!!