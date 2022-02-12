%% Calculo de los desplazamientos y las reacciones en los empotramiento, 
%% las deformaciones y los esfuerzos de una estructura en TENSION PLANA
% NOTA: este codigo SOLO es apropiado para TENSION PLANA usando EFs
% rectangulares serendipitos de 8 nodos

clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% defino las variables/constantes
X = 1; Y = 2; % un par de constantes que ayudaran en la lectura del codigo

%% se define la estructura a calcular
%filename = {'malla_1', 'malla1'};
%filename = {'malla_2', 'malla2'};
%filename = {'malla_3', 'malla3'};
filename = {'malla_4', 'malla4'};
archivo_xlsx = fullfile('..', filename{1}, [filename{2} '.xlsx']);

%% se leen las coordenadas de los nodos
T = leer_excel(archivo_xlsx, 'xnod');
idxNODO         = T{:,'nodo'};
xnod(idxNODO,:) = T{:,{'x','y'}}; % = [x, y]
nno             = size(xnod,1);   % numero de nodos

%% se lee la matriz de conectividad (LaG) y el tipo de material del EF e
T = leer_excel(archivo_xlsx, 'LaG_mat');
idxEF        = T{:,'EF'};
LaG(idxEF,:) = T{:,{'NL1','NL2','NL3','NL4','NL5','NL6','NL7','NL8'}};
mat          = T{:,'material'};% tipo de material para cada EF
nef          = size(LaG,1);    % numero de EFs (numero de filas de LaG)

%% Se definen los grados de libertad
ngdl = 2*nno;        % numero de grados de libertad (dos por nodo)
gdl  = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad

%% definición de los materiales
T = leer_excel(archivo_xlsx, 'prop_mat');
E    = T{:,'E'};       % modulo de elasticidad del solido
nu   = T{:,'nu'};      % coeficiente de Poisson
rho  = T{:,'rho'};     % densidad
t    = T{:,'espesor'}; % espesor
nmat = size(E,1);      % numero de materiales

%% se definen los apoyos y sus desplazamientos
T = leer_excel(archivo_xlsx, 'restric');
idxNODO = T{:,'nodo'};
dirdesp = T{:,'direccion'};
ac      = T{:,'desplazamiento'}; % desplazamientos conocidos en los apoyos

%% Se definen las restricciones 
%{
ngdl_res = size(ac,1); % numero de grados de libertad restringidos
restric = zeros(ngdl_res,2);
for i = 1:ngdl_res
%                       nodo        direccion   desplazamiento    
   restric(i,:) = [ gdl(idxNODO(i), dirdesp(i)) ac(i) ];
end
%}
restric = [ gdl(sub2ind([nno 2], idxNODO, dirdesp)) ac ];

%% se definen las cargas puntuales
T = leer_excel(archivo_xlsx, 'carga_punt');
idxNODO = T{:,'nodo'};
dirfuer = T{:,'direccion'};
f_punt  = T{:,'fuerza_puntual'}; % desplazamientos conocidos en los apoyos

f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
%{
nf_punt = size(f_punt,1); % numero de fuerzas puntuales
for i = 1:nf_punt
%        nodo        direccion      fuerza puntual
   f(gdl(idxNODO(i), dirfuer(i))) = f_punt(i);
end
%}
f(gdl(sub2ind([nno 2], idxNODO, dirfuer))) = f_punt;

%% se leen algunas variables
T        = readcell(archivo_xlsx, 'Sheet','varios','Range','B1:B9');
%E        = T{1}; % modulo de elasticidad E
%nu       = T{2}; % coeficiente de Poisson
%rho      = T{3}; % densidad del material
g        = T{4}; % aceleracion de la gravedad
%t        = T{5}; % espesor de la estructura
U_LONG   = T{6}; % unidades de longitud
U_FUERZA = T{7}; % unidades de fuerza
U_ESFUER = T{8}; % unidades de esfuerzo
ESC_UV   = T{9}; % factor de escala para los desplazamientos

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(LaG(e,[1 3 5 7]),X));
   cgy(e) = mean(xnod(LaG(e,[1 3 5 7]),Y));
   text(cgx(e), cgy(e), num2str(e), 'Color', 'r');
end
plot(xnod(:,X), xnod(:,Y), 'ro');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');

%% Funciones de forma serendipitas del elemento rectangular de 8 nodos:
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa deduccion_funciones_forma/FF_serendipitos_Q4_Q8.m

Nforma = @(xi,eta) [ ...
-((eta - 1)*(xi - 1)*(eta + xi + 1))/4       % N1
((xi^2 - 1)*(eta - 1))/2                     % N2
((eta - 1)*(xi + 1)*(eta - xi + 1))/4        % N3
-((eta^2 - 1)*(xi + 1))/2                    % N4
((eta + 1)*(xi + 1)*(eta + xi - 1))/4        % N5
-((xi^2 - 1)*(eta + 1))/2                    % N6
((eta + 1)*(xi - 1)*(xi - eta + 1))/4        % N7
((eta^2 - 1)*(xi - 1))/2                 ];  % N8

%% Derivadas de N con respecto a xi
dN_dxi = @(xi,eta) [ ...
-((eta + 2*xi)*(eta - 1))/4                  % dN1_dxi
eta*xi - xi                                  % dN2_dxi
((eta - 2*xi)*(eta - 1))/4                   % dN3_dxi
1/2 - eta^2/2                                % dN4_dxi
((eta + 2*xi)*(eta + 1))/4                   % dN5_dxi
-xi*(eta + 1)                                % dN6_dxi
-((eta - 2*xi)*(eta + 1))/4                  % dN7_dxi
eta^2/2 - 1/2                            ];  % dN8_dxi

%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ...
-((2*eta + xi)*(xi - 1))/4                   % dN1_deta
xi^2/2 - 1/2                                 % dN2_deta
((xi + 1)*(2*eta - xi))/4                    % dN3_deta
-eta*(xi + 1)                                % dN4_deta
((2*eta + xi)*(xi + 1))/4                    % dN5_deta
1/2 - xi^2/2                                 % dN6_deta
-((xi - 1)*(2*eta - xi))/4                   % dN7_deta
eta*(xi - 1)                             ];  % dN8_deta

%% Parametros de la cuadratura de Gauss-Legendre
% se calculan las raices x_gl y los pesos w_gl de polinomios de Legendre
n_gl         = 2;                        % orden de la cuadratura
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl);   % matriz de rigidez global como RALA (sparse)
N = cell(nef,n_gl,n_gl); % contenedor para las matrices de forma
B = cell(nef,n_gl,n_gl); % contenedor para las matrices de deformacion

De = cell(nmat,1);
be = cell(nmat,1);
for i = 1:nmat
    % matriz constitutiva del elemento para TENSION PLANA
    De{i} = (E(i)/(1-nu(i)^2)) * [ 1     nu(i)  0
                                  nu(i)  1      0
                                  0      0      (1-nu(i))/2 ];
    be{i} = [0; -rho(i)*g];  % vector de fuerzas masicas del elemento
end

idx = cell(nef,1);
for e = 1:nef          % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:) ...
              gdl(LaG(e,3),:)  gdl(LaG(e,4),:) ...
              gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8),:) ];
   
   % Calculo las matrices de rigidez y el vector de fuerzas nodales
   % equivalentes del elemento
   Ke = zeros(16);
   fe = zeros(16,1);
   det_Je = zeros(n_gl,n_gl); % en esta matriz se almacenaran los Jacobianos

   % Se determinan las coordenadas de los nodos el EF e
   xe = xnod(LaG(e,:),X);
   ye = xnod(LaG(e,:),Y);

   for p = 1:n_gl
      for q = 1:n_gl
         xi_gl  = x_gl(p);
         eta_gl = x_gl(q);
         
         % Se evaluan las funciones de forma en los puntos de integracion
         % de Gauss-Legendre
         NNforma = Nforma(xi_gl, eta_gl);
         
         % Se evaluan las derivadas de las funciones de forma en los puntos
         % de integracion de Gauss-Legendre
         ddN_dxi  = dN_dxi (xi_gl, eta_gl);
         ddN_deta = dN_deta(xi_gl, eta_gl);
         
         dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
         dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
         
         % Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];
            
         % Se calcula el determinante del Jacobiano
         det_Je(p,q) = det(Je);
         
         N{e,p,q} = zeros(2,2*8);
         B{e,p,q} = zeros(3,2*8);
         for i = 1:8
            % Se ensambla la matriz de funciones de forma N
            N{e,p,q}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                        0           NNforma(i) ];
         
            % Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je(p,q);
            dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je(p,q);
            B{e,p,q}(:,[2*i-1 2*i]) = [ dNi_dx 0          % aqui se ensambla
                                        0      dNi_dy     % y asigna la matriz
                                        dNi_dy dNi_dx ];  % B_i
         end

         % se arma la matriz de rigidez del elemento e
         Ke = Ke + B{e,p,q}'*De{mat(e)}*B{e,p,q}*det_Je(p,q)*t(mat(e))*w_gl(p)*w_gl(q);

         % vector de fuerzas nodales equivalentes asociado a la fuerza
         % masica
         fe = fe + N{e,p,q}'*be{mat(e)}*det_Je(p,q)*t(mat(e))*w_gl(p)*w_gl(q);
      end
   end
   
   if any(any(det_Je <= 0))
      error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
   end

   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)      + fe;
end

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero');

%% se leen las cargas distribuidas
T       = leer_excel(archivo_xlsx, 'carga_distr');
idxELEM = T{:,'elemento'};
nodoijk = T{:,{'nodo_i','nodo_j','nodo_k'}};
carga   = T{:,{'tix','tiy','tjx','tjy','tkx','tky'}};
nlcd    = size(carga,1); % numero de lados con carga distribuida

%% Relacion de las cargas superficiales (vector ft)
ft = sparse(ngdl,1); % fuerzas nodales equivalentes de cargas superficiales
for i = 1:nlcd
   e   = idxELEM(i);
   LaG_e = LaG(e,1:8);
   fte = t2ft_Q8(xnod(LaG_e,[X Y]), LaG_e, nodoijk(i,:), carga(i,:), t(mat(e)));
   ft(idx{e},:) = ft(idx{e},:) + fte;
end

% Agrego al vector de fuerzas nodales equivalentes las fuerzas
% superficiales calculadas
f = f + ft;

%% se leen las constantes de balastro k (cimentacion elastica)
T       = leer_excel(archivo_xlsx, 'kWinkler');
idxELEM = T{:,'elemento'};
nodoijk = T{:,{'nodo_i','nodo_j','nodo_k'}};
kwinkl  = T{:,{'kix','kiy','kjx','kjy','kkx','kky'}};
nlk     = size(carga,1); % numero de lados con cimentacion elastica

%% Relacion de los EFs sobre cimentacion elastica
for i = 1:length(idxELEM)
   e = idxELEM(i);
   LaG_e = LaG(e,1:8);
   He = Hwinkler_8(xnod(LaG_e,[X Y]), LaG_e, nodoijk(i,:), kwinkl(i,:), t(mat(e)));
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + He;
end

%% grados de libertad del desplazamiento conocidos y desconocidos  
c = restric(:,1);   d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = restric(:,2);   % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);      % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,2,nno)';
xdef = xnod + ESC_UV*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y), 'Color','b'); % original
   line(xdef(LaG(e,[1:8 1]),X), xdef(LaG(e,[1:8 1]),Y), 'Color','r'); % deformada
end
xlabel(['Eje X [' U_LONG ']']);
ylabel(['Eje Y [' U_LONG ']']);
axis equal tight;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside');
title(sprintf('Deformada escalada %d veces', ESC_UV));

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl,n_gl);
esf = cell(nef,n_gl,n_gl);
for e = 1:nef
   ae = a(idx{e});            % desplazamientos de los gdl del elemento e
   
   for pp = 1:n_gl
      for qq = 1:n_gl
         def{e,pp,qq} = B{e,pp,qq}*ae;           % calculo las deformaciones
         esf{e,pp,qq} = De{mat(e)}*def{e,pp,qq}; % calculo los esfuerzos
      end
   end
end

%% Se extrapolan los esfuerzos y las deformaciones a los nodos y se alisan
% adicionalmente se calcula el error en el alisado
[sx,  error_sx ] = extrapolar_esf_def(xnod, LaG, esf, 'sx');
[sy,  error_sy ] = extrapolar_esf_def(xnod, LaG, esf, 'sy');
[txy, error_txy] = extrapolar_esf_def(xnod, LaG, esf, 'txy');
[ex,  error_ex ] = extrapolar_esf_def(xnod, LaG, def, 'ex');
[ey,  error_ey ] = extrapolar_esf_def(xnod, LaG, def, 'ey');
[gxy, error_gxy] = extrapolar_esf_def(xnod, LaG, def, 'gxy');

%% en tension plana ...
sz   = zeros(nno,1);
txz  = zeros(nno,1);
tyz  = zeros(nno,1);

ez   = -(nu/E)*(sx+sy);              % deformaciones ez
tmax = sqrt(((sx-sy)/2).^2+txy.^2);  % esfuerzo cortante maximo
s1   = (sx+sy)/2 + tmax;             % esfuerzo normal maximo
s2   = (sx+sy)/2 - tmax;             % esfuerzo normal minimo
ang  = 0.5*atan2(2*txy, sx-sy);      % angulo de inclinacion de s1

s3   = zeros(nno,1);
sv   = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2); % von Mises

%% Se reportan los resultados en un archivo .xlsx
tabla_aq = array2table([(1:nno)', reshape(a,2,nno)', reshape(q,2,nno)'], ...
    'VariableNames', {'nodo', ['u_' U_LONG], 'v', ['qx_' U_FUERZA], 'qy'});

tabla_def = array2table([(1:nno)', ex, ey, ez, gxy], ...
    'VariableNames', {'nodo', 'ex', 'ey', 'ez', 'gxy_rad'});

tabla_esf = array2table([(1:nno)', sx, sy, txy, s1, s2, ang, tmax, sv], ...
    'VariableNames', {'nodo', ['sx_' U_ESFUER], 'sy', 'txy', 's1', 's2', 'ang_rad', 'tmax', 'sv'});

filename_results = ['resultados_' filename{2} '.xlsx'];
writetable(tabla_aq,  filename_results, 'Sheet', 'aq')
writetable(tabla_def, filename_results, 'Sheet', 'deformaciones')
writetable(tabla_esf, filename_results, 'Sheet', 'esfuerzos')

fprintf('Calculo finalizado. Resultados en "%s".\n', filename_results);

%% se grafican las deformaciones
figure
subplot(1,4,1); plot_def_esf(xnod, LaG, ex,  '\epsilon_x')
subplot(1,4,2); plot_def_esf(xnod, LaG, ey,  '\epsilon_y')
subplot(1,4,3); plot_def_esf(xnod, LaG, ez,  '\epsilon_z')
subplot(1,4,4); plot_def_esf(xnod, LaG, gxy, '\gamma_{xy} [rad]')

%% se grafican los esfuerzos
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, sx,  ['\sigma_x [' U_ESFUER ']'])
subplot(1,3,2); plot_def_esf(xnod, LaG, sy,  ['\sigma_y [' U_ESFUER ']'])
subplot(1,3,3); plot_def_esf(xnod, LaG, txy, ['\tau_{xy} [' U_ESFUER ']'])

%% se grafican los errores en los esfuerzos
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, error_sx,  ['Error \sigma_x [' U_ESFUER ']'])
subplot(1,3,2); plot_def_esf(xnod, LaG, error_sy,  ['Error \sigma_y [' U_ESFUER ']'])
subplot(1,3,3); plot_def_esf(xnod, LaG, error_txy, ['Error \tau_{xy} [' U_ESFUER ']'])

%% se grafican los esfuerzos principales y el esfuerzo cortante maximo
figure
subplot(1,3,1); plot_def_esf(xnod, LaG, s1,   ['(\sigma_1)_{xy} [' U_ESFUER ']'], { ang })
subplot(1,3,2); plot_def_esf(xnod, LaG, s2,   ['(\sigma_2)_{xy} [' U_ESFUER ']'], { ang+pi/2 })
subplot(1,3,3); plot_def_esf(xnod, LaG, tmax, ['\tau_{max} [' U_ESFUER ']'],      { ang+pi/4, ang-pi/4 })

%% se grafican los esfuerzos de von Mises
figure
plot_def_esf(xnod, LaG, sv, ['Esfuerzos de von Mises [' U_ESFUER ']']);

%% se exportan los resultados a GiD/Paraview
% Pasando los esfuerzos ya promediados:
%export_to_GiD('c5_ejemplo_a',xnod,LaG,a,q,[sx sy sz txy txz tyz]);

% Pasando los puntos de Gauss [RECOMENDADO] !!!
% export_to_GiD('c5_ejemplo_b',xnod,LaG,a,q,esf);                    

%%
return; % bye, bye!

%% Lee del archivo "nombre_archivo" de EXCEL la hoja "hoja"
function H = leer_excel(archivo_xlsx, hoja)

    if verLessThan('matlab', '9.9') % R2019b or older
        H = readtable(archivo_xlsx, 'Sheet', hoja);
    else
        H = readtable(archivo_xlsx, 'Sheet', hoja, 'format', 'auto');
    end
end

%% Grafica los esfuerzos y las deformaciones
function plot_def_esf(xnod, LaG, variable, texto, angulos)
    X = 1; Y = 2;
    hold on; 
    colorbar;
    
    nef = size(LaG, 1);    
    for e = 1:nef
       data.numEF = e; 
       data.LaG_e = LaG(e,:);
       fill(xnod(LaG(e,:),X), xnod(LaG(e,:),Y), variable(LaG(e,:)), ...
           'UserData', data);
    end
    axis equal tight
    colormap jet
    title(texto);

    esc = 0.5;
    if nargin == 5
        norma = 1; % = variable % si se quiere proporcional
        for i = 1:length(angulos)
            % se indica la flecha de la direccion principal
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}), norma.*sin(angulos{i}),... 
                esc, ...                  % con una escala esc
                'k',...                   % de color negro
                'ShowArrowHead','off',... % una flecha sin cabeza
                'LineWidth',2, ...        % con un ancho de linea 2
                'Marker','.');            % y en el punto (x,y) poner un punto '.'
            
            % la misma flecha girada 180 grados
            quiver(xnod(:,X),xnod(:,Y),...             
                norma.*cos(angulos{i}+pi), norma.*sin(angulos{i}+pi),... 
                esc,'k', 'ShowArrowHead','off', 'LineWidth',2, 'Marker','.');                    
        end      
    end

    % para mostrar el tooltip
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.UpdateFcn = @mostrar_info_nodo;    
end

function txt = mostrar_info_nodo(~,info)
    if strcmp(info.Target.Type, 'patch')
        x      = info.Position(1);
        y      = info.Position(2);
        numEF  = info.Target.UserData.numEF;
        xnod_e = info.Target.Vertices;
    
        % se busca el punto más cercano
        [~, idx_e] = min(hypot(xnod_e(:,1) - x, xnod_e(:,2) - y));
        numNOD  = info.Target.UserData.LaG_e(idx_e);
        val_esf = info.Target.CData(idx_e);

        txt = ['(x,y) = (' num2str(x) ',' num2str(y) '), ' ...
               'nodo = '     num2str(numNOD)          ', ' ...
               'EF = '       num2str(numEF)           ', ' ...
               'variable = ' num2str(val_esf)];
    else
        txt = 'Haga zoom y evite seleccionar el quiver()';
    end
end

%% Extrapola/alisa esfuerzos y deformaciones de puntos de Gauss a los nodos
function [esf, error_esf] = extrapolar_esf_def(xnod, LaG, esfuerzo, tipo_esf)
    nno = size(xnod, 1);
    nef = size(LaG, 1);

    num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
    esf.sum      = zeros(nno,1);
    esf.max      =  -inf(nno,1);
    esf.min      =   inf(nno,1);

    A = [ ... 
      3^(1/2)/2 + 1,            -1/2,            -1/2,   1 - 3^(1/2)/2
    3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4
               -1/2,   1 - 3^(1/2)/2,   3^(1/2)/2 + 1,            -1/2
    1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4
      1 - 3^(1/2)/2,            -1/2,            -1/2,   3^(1/2)/2 + 1
    1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 3^(1/2)/4 + 1/4
               -1/2,   3^(1/2)/2 + 1,   1 - 3^(1/2)/2,            -1/2
    3^(1/2)/4 + 1/4, 3^(1/2)/4 + 1/4, 1/4 - 3^(1/2)/4, 1/4 - 3^(1/2)/4 ];

    switch tipo_esf
        case {'sx',  'ex'},  num_esf = 1;
        case {'sy',  'ey'},  num_esf = 2;
        case {'txy', 'gxy'}, num_esf = 3;
        otherwise,           error('Opcion no soportada');
    end

    for e = 1:nef
        esf_EF_e = A * [ esfuerzo{e,1,1}(num_esf)
                         esfuerzo{e,1,2}(num_esf)
                         esfuerzo{e,2,1}(num_esf)
                         esfuerzo{e,2,2}(num_esf) ];        
        
        esf.sum(LaG(e,:),:) = esf.sum(LaG(e,:),:)    + esf_EF_e;
        esf.max(LaG(e,:),:) = max(esf.max(LaG(e,:),:), esf_EF_e);
        esf.min(LaG(e,:),:) = min(esf.max(LaG(e,:),:), esf_EF_e);     
                                                
        num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
    end

    %% alisado (promedio de los esfuerzos en los nodos)
    esf.prom = esf.sum./num_elem_ady;    
    
    %% variables a retornar
    error_esf = (esf.max - esf.min)./esf.prom; % error en el alisado
    error_esf = log10(abs(error_esf));
    error_esf(error_esf < log10(0.1)) = -3;
    esf       = esf.prom;                      % esfuerzo promedio
end
