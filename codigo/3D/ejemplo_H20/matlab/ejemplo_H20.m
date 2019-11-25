% clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% -----------------------------------------------------------------------------
%% NOTA: este codigo SOLO es apropiado para EFs hexaedricos de 20 nodos (H20)
%% -----------------------------------------------------------------------------

%% defino las variables/constantes
X    = 1;                % constantes que ayudaran en la lectura del codigo
Y    = 2;
Z    = 3;
Ee   = 200e9;            % modulo de elasticidad del solido (Pa) = 200 GPa
nue  = 0.30;             % coeficiente de Poisson
rhoe = 7850;             % densidad (kg/m^3)
g    = 9.81;             % aceleracion de la gravedad (m/s^2)
be   = [0; -rhoe*g; 0];  % vector de fuerzas masicas del elemento (el eje Y es el vertical)

%% cargar
% xnod - posicion de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
%filename = 'mallas/malla_H20_viga';
filename = 'mallas/malla_H20_conexion';
[xnod, LaG]            = import_from_GiD(filename);
[cargas,restricciones] = import_cargas_restric(filename);

nno  = size(xnod,1);     % numero de nodos (numero de filas de xnod)
nef  = size(LaG,1);      % numero de EFs (numero de filas de LaG)
nnpe = size(LaG,2);      % numero de nodos por elemento finito (= 20)
ngdl = 3*nno;            % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

if nnpe ~= 20
   error('Este codigo SOLO es sirve para EFs hexaedricos de 20 nodos (H20)')
end

%% Se definen las restricciones 
ngdl_res = size(restricciones,1); % numero de grados de libertad restringidos
restric = zeros(ngdl_res,2);
for i = 1:ngdl_res
%                       nodo                direccion           desplazamiento    
   restric(i,:) = [ gdl(restricciones(i,1), restricciones(i,2)) restricciones(i,3) ];
end

%% Relacion de cargas puntuales
f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
ncargas = size(cargas,1);
for i = 1:ncargas
%        nodo         direccion       desplazamiento    
   f(gdl(cargas(i,1), cargas(i,2))) = cargas(i,3);
end

%% Se dibuja la malla de elementos finitos
figure;
hold on;
for e = 1:nef
   % Se está graficando del mismo modo que GiD lo hace Y es el vertical
   plot_EF_H20(e, xnod, LaG);
end
daspect([1 1 1]);
view(3)
title('Malla de elementos finitos','FontSize',26);
xlabel('Eje Z');
ylabel('Eje X');
zlabel('Eje Y');
grid on;
box on;

%% Funciones de forma del elemento H20:
Nforma   = @N_H20;
dN_dxi   = @dN_dxi_H20;
dN_deta  = @dN_deta_H20;
dN_dzeta = @dN_dzeta_H20;

%% Parametros de la cuadratura de Gauss-Legendre para hexaedros
[x_gl, w_gl] = gausslegendre_quad_hexa(2); % 2 puntos de integracion por eje
n_gl = length(w_gl);

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl, ngdl);  % matriz de rigidez global como RALA (sparse)
N = cell(nef, n_gl);     % contenedor para las matrices de forma
B = cell(nef, n_gl);     % contenedor para las matrices de deformacion

% matriz constitutiva del elemento tridimensional
d1 = (1-nue)/(1-2*nue);
d2 = nue/(1-2*nue);
De = Ee/(1+nue)*[ ...
   d1 d2 d2 0   0   0
   d2 d1 d2 0   0   0
   d2 d2 d1 0   0   0
   0  0  0  1/2 0   0
   0  0  0  0   1/2 0
   0  0  0  0   0   1/2 ];

for e = 1:nef          % ciclo sobre todos los elementos finitos
   % Calculo la matriz de rigidez y el vector de fuerzas nodales
   % equivalentes del elemento
   Ke = zeros(3*nnpe);
   fe = zeros(3*nnpe, 1);
   det_Je = zeros(n_gl, 1); % en esta matriz se almacenaran los Jacobianos
   
   for p = 1:n_gl
      xi_gl   = x_gl(p,1);
      eta_gl  = x_gl(p,2);
      zeta_gl = x_gl(p,3);
      
      % Se evaluan las funciones de forma en los puntos de integracion
      % de Gauss-Legendre
      NNforma = Nforma(xi_gl, eta_gl, zeta_gl);
      
      % Se evaluan las derivadas de las funciones de forma en los puntos
      % de integracion de Gauss-Legendre
      ddN_dxi   = dN_dxi  (xi_gl, eta_gl, zeta_gl);  xe = xnod(LaG(e,:),X);
      ddN_deta  = dN_deta (xi_gl, eta_gl, zeta_gl);  ye = xnod(LaG(e,:),Y);
      ddN_dzeta = dN_dzeta(xi_gl, eta_gl, zeta_gl);  ze = xnod(LaG(e,:),Z);
      
      dx_dxi   = sum(ddN_dxi  .*xe);   dy_dxi   = sum(ddN_dxi  .*ye);   dz_dxi   = sum(ddN_dxi  .*ze);
      dx_deta  = sum(ddN_deta .*xe);   dy_deta  = sum(ddN_deta .*ye);   dz_deta  = sum(ddN_deta .*ze);
      dx_dzeta = sum(ddN_dzeta.*xe);   dy_dzeta = sum(ddN_dzeta.*ye);   dz_dzeta = sum(ddN_dzeta.*ze);
      
      % Se ensambla la matriz Jacobiana del elemento
      Je = [ dx_dxi    dy_dxi    dz_dxi
             dx_deta   dy_deta   dz_deta
             dx_dzeta  dy_dzeta  dz_dzeta ];
      
      % Se calcula el determinante del Jacobiano
      det_Je(p) = det(Je);
      
      N{e,p} = zeros(3,3*nnpe);
      B{e,p} = zeros(6,3*nnpe);
      for i = 1:nnpe
         % Se ensambla la matriz de funciones de forma N
         N{e,p}(:,3*i-[2 1 0]) = [ ...
            NNforma(i)  0           0
            0           NNforma(i)  0
            0           0           NNforma(i) ];
         
         % Se ensambla la matriz de deformacion del elemento B
         tmp = Je\[ ddN_dxi(i);  ddN_deta(i);  ddN_dzeta(i) ];
         dNi_dx = tmp(1);
         dNi_dy = tmp(2);
         dNi_dz = tmp(3);
         B{e,p}(:, 3*i-[2 1 0]) = [ ...
            dNi_dx 0      0        % aqui se ensambla
            0      dNi_dy 0        % y asigna la matriz
            0      0      dNi_dz   % B_i
            dNi_dy dNi_dx 0
            dNi_dz 0      dNi_dx
            0      dNi_dz dNi_dy ];
      end
      
      % se arma la matriz de rigidez del elemento e
      Ke = Ke + B{e,p}'*De*B{e,p}*det_Je(p)*w_gl(p);
      
      % vector de fuerzas nodales equivalentes
      fe = fe + N{e,p}'*be*det_Je(p)*w_gl(p);
   end
   
   if any(det_Je <= 0)
      error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
   end
   
   % Se sacan los grados de libertad del hexaedro
   idx = zeros(1,3*nnpe);
   for i = 1:nnpe
      idx(3*i-[2 1 0]) = gdl(LaG(e,i),:);
   end
   
   K(idx,idx) = K(idx,idx) + Ke;
   f(idx,:)   = f(idx,:)   + fe;
   % NOTA: para quitar el Warning de MATLAB que el indexing is slow mire:
   % https://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/
end

%% Muestro la configuracion de la matriz K (K es rala)
figure;
spy(K);
title('Los puntos representan los elementos diferentes de cero', 'FontSize', 26);

%% grados de libertad del desplazamiento conocidos y desconocidos
c = restric(:,1);
d = setdiff(1:ngdl,c)';

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
ac = restric(:,2); % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% imprimo los resultados
%format short g
%disp('Nodo   Despl_x (m)   Despl_y (m)   Despl_z (m)= ');  [1:nno; reshape(a,3,nno)]'
%disp('Nodo Fuerzas nodales equiv. X, Y, Z (N) = ');        [1:nno; reshape(f,3,nno)]'
%disp('Nodo Fuerzas nodales equil. X, Y, Z (N) = ');        [1:nno; reshape(q,3,nno)]'

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,3,nno)';
escala = 1000;             % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   % Se está graficando del mismo modo que GiD lo hace Y es el vertical
   
   % original
   h1 = plot_EF_H20(e, xnod, LaG);
   set(h1, 'Color', [0 0 1]); % color expresado en notacion RBG entre 0 y 1
   
   % deformada
   h2 = plot_EF_H20(e, xdef, LaG);
   set(h2, 'Color', [1 0 0]);
end
daspect([1 1 1]);
view(3)
title('Malla de elementos finitos','FontSize',26);
xlabel('Eje Z');
ylabel('Eje X');
zlabel('Eje Y');
grid on;
box on;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside')
title(sprintf('Deformada escalada %d veces',escala), 'FontSize', 26);

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl);
esf = cell(nef,n_gl);
for e = 1:nef
   idx = zeros(1,3*nnpe);
   for i = 1:nnpe
      idx(3*i-[2 1 0]) = gdl(LaG(e,i),:);
   end
   
   ae = a(idx);                % desplazamientos de los gdl del elemento e
   
   for p = 1:n_gl
      def{e,p} = B{e,p}*ae;    % calculo las deformaciones
      esf{e,p} = De*def{e,p};  % calculo los esfuerzos
   end
end

%% Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
sx  = zeros(nno,1);       ex  = zeros(nno,1);
sy  = zeros(nno,1);       ey  = zeros(nno,1);
sz  = zeros(nno,1);       ez  = zeros(nno,1);
txy = zeros(nno,1);       gxy = zeros(nno,1);
txz = zeros(nno,1);       gxz = zeros(nno,1);
tyz = zeros(nno,1);       gyz = zeros(nno,1);

A = matriz_extrapolacion_esfuerzos_H20();

v_esf = zeros(n_gl,1);
v_def = zeros(n_gl,1);
for e = 1:nef
   for p = 1:n_gl, v_esf(p) = esf{e,p}(1); end;   sx(LaG(e,:), :) =  sx(LaG(e,:),:) + A*v_esf;
   for p = 1:n_gl, v_esf(p) = esf{e,p}(2); end;   sy(LaG(e,:), :) =  sy(LaG(e,:),:) + A*v_esf;
   for p = 1:n_gl, v_esf(p) = esf{e,p}(3); end;   sz(LaG(e,:), :) =  sz(LaG(e,:),:) + A*v_esf;
   for p = 1:n_gl, v_esf(p) = esf{e,p}(4); end;   txy(LaG(e,:),:) = txy(LaG(e,:),:) + A*v_esf;
   for p = 1:n_gl, v_esf(p) = esf{e,p}(5); end;   txz(LaG(e,:),:) = txz(LaG(e,:),:) + A*v_esf;
   for p = 1:n_gl, v_esf(p) = esf{e,p}(6); end;   tyz(LaG(e,:),:) = tyz(LaG(e,:),:) + A*v_esf;
   
   for p = 1:n_gl, v_def(p) = def{e,p}(1); end;   ex(LaG(e,:), :) =  ex(LaG(e,:),:) + A*v_def;
   for p = 1:n_gl, v_def(p) = def{e,p}(2); end;   ey(LaG(e,:), :) =  ey(LaG(e,:),:) + A*v_def;
   for p = 1:n_gl, v_def(p) = def{e,p}(3); end;   ez(LaG(e,:), :) =  ez(LaG(e,:),:) + A*v_def;
   for p = 1:n_gl, v_def(p) = def{e,p}(4); end;   gxy(LaG(e,:),:) = gxy(LaG(e,:),:) + A*v_def;
   for p = 1:n_gl, v_def(p) = def{e,p}(5); end;   gxz(LaG(e,:),:) = gxz(LaG(e,:),:) + A*v_def;
   for p = 1:n_gl, v_def(p) = def{e,p}(6); end;   gyz(LaG(e,:),:) = gyz(LaG(e,:),:) + A*v_def;
   
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;        ex  =  ex./num_elem_ady;
sy  =  sy./num_elem_ady;        ey  =  ey./num_elem_ady;
sz  =  sz./num_elem_ady;        ez  =  ez./num_elem_ady;
txy = txy./num_elem_ady;        gxy = gxy./num_elem_ady;
txz = txz./num_elem_ady;        gxz = gxz./num_elem_ady;
tyz = tyz./num_elem_ady;        gyz = gyz./num_elem_ady;

%% Se imprimen las deformaciones en los nodos
%disp('Deformaciones: (Nodo,ex,ey,ez,gxy,gxz,gyz) = ');
%disp([(1:nno)'  ex  ey  ez  gxy  gxz  gyz])

%% Se imprimen los esfuerzos en los nodos
%disp('Esfuerzos (Pa):  (Nodo,sx,sy,sz,txy,txz,tyz) = ');
%disp([(1:nno)'  sx  sy  sz  txy  txz  tyz])

%% Se calculan para cada nodo los esfuerzos principales y sus direcciones
s1 = zeros(nno,1);  dir_s1 = zeros(nno, 3);
s2 = zeros(nno,1);  dir_s2 = zeros(nno, 3);
s3 = zeros(nno,1);  dir_s3 = zeros(nno, 3);
for i = 1:nno
   [dirppales, esfppales] = eig([sx(i)   txy(i)  txz(i)    % matriz de esfuerzos
                                 txy(i)  sy(i)   tyz(i)    % de Cauchy
                                 txz(i)  tyz(i)  sz(i)]);

   [esfppales,idx] = sort(diag(esfppales), 'descend');  
   s1(i) = esfppales(1);   dir_s1(i,:) = dirppales(idx(1),:)';
   s2(i) = esfppales(2);   dir_s2(i,:) = dirppales(idx(2),:)';
   s3(i) = esfppales(3);   dir_s3(i,:) = dirppales(idx(3),:)';
end

% Esfuerzo cortante maximo
tmax = (s1-s3)/2;                               % esfuerzo cortante maximo
   
%% Calculo de los esfuerzos de von Mises
sv = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2);

%% imprimo los resultados
%disp('Nodo,s1(Pa),s2(Pa), s3(Pa), tmax(Pa), Esfuerzos de von Mises (Pa) = ');
%disp([(1:nno)'  s1  s2  s3  tmax  sv])

% Pasando los esfuerzos ya promediados:
export_to_GiD('resultados/conexion_esf_nodos',xnod,LaG,a,q,[sx sy sz txy txz tyz]);

%% Pasando los puntos de Gauss [RECOMENDADO] !!!
export_to_GiD('resultados/conexion_H20_esf_GP',xnod,LaG,a,q,esf);


%% Se genera un archivo .VTK para visualizar en Paraview
%{
meshio.write_points_cells(
    "resultados.vtk",
    points=xnod,
    cells={"quad8": LaG[:,[0,2,4,6,1,3,5,7]] },
    point_data = {
        'ex':ex, 'ey':ey, 'ez':ez,     'gxy':gxy,
        'sx':sx, 'sy':sy, 'txy':txy,
        's1':s1, 's2':s2, 'tmax':tmax, 'sv':sv,
        'uv'  :a.reshape((nno,2)),
        's1n1':np.c_[s1*np.cos(ang),           s1*np.sin(ang)          ],
        's2n2':np.c_[s2*np.cos(ang + np.pi/2), s2*np.sin(ang + np.pi/2)]
        }
    # cell_data=cell_data,
    # field_data=field_data
)
%}

%{
vtkwrite('resultados/resultados.vtk', 'structured_grid', ...
    xnod(:,X), xnod(:,Y), xnod(:,Z), ...
    'scalars', 's1',  s1,  ...
    'scalars', 's2',  s2,  ...
    'scalars', 's3',  s3,  ...        
    'scalars', 'sv',  sv,  ...            
    'scalars', 'sx',  sx,  ...
    'scalars', 'sy',  sy,  ...
    'scalars', 'sz',  sz,  ...    
    'scalars', 'txy', txy, ...
    'scalars', 'txz', txz, ...    
    'scalars', 'tyz', tyz, ...        
    'scalars', 'tmax', tmax, ...          
    'scalars', 'ex',  ex,  ...
    'scalars', 'ey',  ey,  ...
    'scalars', 'ez',  ez,  ...    
    'scalars', 'gxy', gxy, ...    
    'scalars', 'gxz', gxz, ...    
    'scalars', 'gyz', gyz, ...
    'vectors','uvw', delta(:,X), delta(:,Y), delta(:,Z), ...
    'vectors','n1',   dir_s1(:,X), dir_s1(:,Y), dir_s1(:,Z), ...    
    'vectors','n2',   dir_s2(:,X), dir_s2(:,Y), dir_s2(:,Z), ...
    'vectors','n3',   dir_s3(:,X), dir_s3(:,Y), dir_s3(:,Z), ...  
    'vectors','s1n1', s1.*dir_s1(:,X), s1.*dir_s1(:,Y), s1.*dir_s1(:,Z), ...    
    'vectors','s2n2', s2.*dir_s2(:,X), s2.*dir_s2(:,Y), s2.*dir_s2(:,Z), ...
    'vectors','s3n3', s3.*dir_s3(:,X), s3.*dir_s3(:,Y), s3.*dir_s3(:,Z));
%}

%%
disp('Fin del calculo')
return; % bye, bye!
