function [Mxast_sup, Myast_sup, Mxast_inf, Myast_inf] = ...
                              refuerzo_losa_WoodArmer(Mx, My, Mxy)
% Calculo los momentos de flexión Mxast_sup y Myast_sup 
% asociados al refuerzo en la parte superior de la losa:
if Mx >= -abs(Mxy) && My >= -abs(Mxy) % Caso 1
    Mxast_sup = Mx + abs(Mxy);
    Myast_sup = My + abs(Mxy);
elseif Mx*My > Mxy^2                  % Caso 7
    Mxast_sup = 0;
    Myast_sup = 0;
elseif Mx < -abs(Mxy)                 % Caso 3
    Mxast_sup = 0;
    Myast_sup = My - Mxy^2/Mx;
elseif My < -abs(Mxy)                 % Caso 4
    Mxast_sup = Mx - Mxy^2/My;
    Myast_sup = 0;
else
    error('Error en el cálculo')
end 

% Calculo los momentos de flexión Mxast_inf y Myast_inf 
% asociados al refuerzo en la parte inferior de la losa:
if Mx <= abs(Mxy) && My <= abs(Mxy)   % Caso 2
    Mxast_inf = Mx - abs(Mxy);
    Myast_inf = My - abs(Mxy);
elseif Mx*My > Mxy^2                  % Caso 8
    Mxast_inf = 0;     
    Myast_inf = 0;
elseif Mx > abs(Mxy)                  % Caso 5
    Mxast_inf = 0;
    Myast_inf = My - Mxy^2/Mx;
elseif My > abs(Mxy)                  % Caso 6
    Mxast_inf = Mx - Mxy^2/My;
    Myast_inf = 0;
else
    error('Error en el cálculo')
end
