function feloc = feloc_EE(w, x1, y1, x2, y2)
% Calcula el vector de fuerzas nodales equivalentes en coord locales para
% una barra sometida a peso propio w
% empotrada a la izquierda y
% empotrada a la derecha

%% Se calcula la longitud de la barra
L = hypot(x2-x1, y2-y1);  

%% Se calculan las cargas distribuidas wx y wy
c = (x2-x1)/L;   s = (y2-y1)/L;  % sin y cos de la inclinacion
wx = w*s;
wy = w*c;

%% se calcula el vector de fuerzas nodales equivalentes
feloc = [wx*L/2;  wy*L/2;  +wy*L^2/12;  wx*L/2;  wy*L/2;  -wy*L^2/12 ];
return