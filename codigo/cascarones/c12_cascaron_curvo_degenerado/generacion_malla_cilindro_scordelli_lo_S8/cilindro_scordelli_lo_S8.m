E  = 4.32e8;         % modulo de elasticidad del solido (Pa)
nu = 0;              % coeficiente de Poisson
qdistr = -90;        % carga (N/m^2)

%% Definimos la geometria de la losa 
% Se genero con generacion_malla_cilindro_scordelli_lo.m
r24_875  = load('scordelli_lo_curvo_S8_24_875', 'LaG','xnod');
r25_125 = load('scordelli_lo_curvo_S8_25_125', 'LaG','xnod');
if ~isequal(r24_875.LaG, r25_125.LaG)
   error('Las matrices LaG no coinciden');
end
LaG = r24_875.LaG;
xnod_lo = r24_875.xnod;
xnod_up = r25_125.xnod;
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