function [Mxast_sup, Myast_sup, Mxast_inf, Myast_inf] = WoodArmer(Mx, My, Mxy)
% Calculo los momentos de flexion Mxast_sup y Myast_sup 
% asociados al refuerzo en la parte superior de la losa:
% Se aplican las ecuaciones /*@\eqref{eq:Mxast_Myast_caso1_WA}@*/
Mxast_sup = Mx + abs(Mxy);
Myast_sup = My + abs(Mxy);
if Mxast_sup < 0 && Myast_sup < 0
    % no se requiere refuerzo en la parte superior de la losa
    Mxast_sup = 0;
    Myast_sup = 0;
else
    if Mxast_sup < 0
        Mxast_sup = 0;
        Myast_sup = My + abs(Mxy^2/Mx);
        if Myast_sup < 0
            Myast_sup = 0;
        end
    end    
    if Myast_sup < 0
        Mxast_sup = Mx + abs(Mxy^2/My);
        Myast_sup = 0;
        if Mxast_sup < 0
            Mxast_sup = 0;
        end        
    end
end 

% Calculo los momentos de flexion Mxast_inf y Myast_inf 
% asociados al refuerzo en la parte inferior de la losa:
% Se aplican las ecuaciones /*@\eqref{eq:Mxast_Myast_caso2_WA}@*/
Mxast_inf = Mx - abs(Mxy);
Myast_inf = My - abs(Mxy);
if Mxast_inf > 0 && Myast_inf > 0
    % no se requiere refuerzo en la parte inferior de la losa
    Mxast_inf = 0;
    Myast_inf = 0;
else
    if Mxast_inf > 0
        Mxast_inf = 0;
        Myast_inf = My - abs(Mxy^2/Mx);
        if Myast_inf > 0
            Myast_inf = 0;
        end        
    end    
    if Myast_inf > 0
        Mxast_inf = Mx - abs(Mxy^2/My);
        Myast_inf = 0;
        if Mxast_inf > 0
            Mxast_inf = 0;
        end                
    end
end

return 