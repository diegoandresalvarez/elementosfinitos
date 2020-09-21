function feloc = feloc_ER(w, x1, y1, x2, y2)
% Calcula el vector de fuerzas nodales equivalentes en coord locales para
% una barra sometida a peso propio w
% empotrada a la izquierda y
% con rotula a la derecha

%% Se calcula la longitud de la barra
L = hypot(x2-x1, y2-y1);

%% Se calculan las cargas distribuidas wx y wy
c = (x2-x1)/L;   s = (y2-y1)/L;  % sin y cos de la inclinacion
wx = w*s;
wy = w*c;

%% se calcula el vector de fuerzas nodales equivalentes
feloc = [wx*L/2;  5*wy*L/8;  wy*L^2/8;  wx*L/2;  3*wy*L/8 ];
return

% Se calcularon con el programa calcular_carga_nodal_equivalente_w_ER.m
% Y1 = (5*L*w)/8 
% Y2 = (3*L*w)/8 
% M1 = (L^2*w)/8 
% M2 = 0
