function out = reinforcement_slab_shell(m11, m22, m12, f11, f22, f12, ...
                                                      fy, h, Ct1, Ct2, Cb1, Cb2)

% NOTA: falta verificar este codigo
%
% Calculation of the slab/shell reinfocement acoording to the sandwich model 
% proposed in the EUROCODE 2 Design of concrete structures: DD ENV 1992-1-1 1992
%
% This is based on:
% * EUROCODE 2 Design of concrete structures: DD ENV 1992-1-1 1992 (sect. A2.8)
% * COMPUTERS AND STRUCTURES INC. Concrete shell reinforcement design. This 
%   document is available in the SAP2000/ETABS documentation.
%
% INPUT: 
% m11, m22, m12 % moments
% f11, f22, f12 % membrane forces: 
% fy            % yield strength
% h             % shell thickness
% Ct1           % thickness of top    cover reinforcement in dir 1 
% Ct2           % thickness of top    cover reinforcement in dir 2
% Cb1           % thickness of bottom cover reinforcement in dir 1 
% Cb2           % thickness of bottom cover reinforcement in dir 2
%
% OUTPUT:
% out is a structure that contains the following fields:
% N11_top, N22_top, N11_bot, N22_bot         % equivalent membrane forces (per unit length)
% N12_top, N12_bot                           % equivalent in-plane forces (per unit length)
% NDes1_top, NDes2_top, NDes1_bot, NDes2_bot % equivalent reinforcement forces (per unit length)
% Ast1_top,  Ast1_bot,  Ast2_top,  Ast2_bot  % rebar area (area per unit length)
% Fc_top, Fc_bot                             % concrete compresive force (per unit length)
% Sc_top, Sc_bot                             % concrete compresive stress (stress units)
%
% By:
% Diego Andrés Alvarez Marín (daalvarez@unal.edu.co)

% follow the convention for the moments in Diego's course:
m11 = -m11;
m22 = -m22;
m12 = -m12;

%% Calculation of distances
dt1 = h/2 - Ct1;                       db1 = h/2 - Cb1;
dt2 = h/2 - Ct2;                       db2 = h/2 - Cb2;
dtmin = min(dt1, dt2);                 dbmin = min(db1, db2);

d1 = h - Ct1 - Cb1;  % = dt1 + db1 = distance between bars in dir 1
d2 = h - Ct2 - Cb2;  % = dt2 + db2 = distance between bars in dir 2
dmin = min(d1, d2);

%% Equivalent membrane forces (per unit length)
% top layer                            % bottom layer
N11_top = (-m11 + f11*db1)/d1;         N11_bot = (+m11 + f11*dt1)/d1;
N22_top = (-m22 + f22*db2)/d2;         N22_bot = (+m22 + f22*dt2)/d2;

%% Equivalent in plane shear forces (per unit length)
% top layer                            % bottom layer
N12_top = (-m12 + f12*dbmin)/dmin;     N12_bot = (+m12 + f12*dtmin)/dmin;

%% Reinforcement in the top and bottom layer
[NDes1_top, NDes2_top, Fc_top] = arrayfun(@calc_NDes_Fc, N11_top, N22_top, N12_top);
[NDes1_bot, NDes2_bot, Fc_bot] = arrayfun(@calc_NDes_Fc, N11_bot, N22_bot, N12_bot);

%% Reinforcement intensities = rebar area per unit width (m^2/m)
phi_s = 0.9;                            % stress reduction factor
Ast1_top = NDes1_top/(phi_s*fy);
Ast2_top = NDes2_top/(phi_s*fy);
Ast1_bot = NDes1_bot/(phi_s*fy);
Ast2_bot = NDes2_bot/(phi_s*fy);

% The thickness of each outer layer is the minumum of
% a) twice the cover measured to the center of the outer reinforcement
thick_t1 = 2*min(Ct1,Ct2); 
thick_b1 = 2*min(Cb1,Cb2); 

% b) twice the distance from the center of the slab to the center of the other reinforcement
thick_t2 = 2*(h/2 - min(Ct1,Ct2));
thick_b2 = 2*(h/2 - min(Cb1,Cb2));

tt = min(thick_t1, thick_t2);
tb = min(thick_b1, thick_b2);

%% Concrete compresive stresses in: [stress units]
Sc_top = Fc_top/tt;  % top layer
Sc_bot = Fc_bot/tb;  % bottom layer

% Outputs:
out.N11_top   = N11_top;
out.N22_top   = N22_top;
out.N11_bot   = N11_bot;
out.N22_bot   = N22_bot;
out.N12_top   = N12_top;
out.N12_bot   = N12_bot;
out.NDes1_top = NDes1_top;
out.NDes2_top = NDes2_top;
out.NDes1_bot = NDes1_bot;
out.NDes2_bot = NDes2_bot;
out.Ast1_top  = Ast1_top;
out.Ast1_bot  = Ast1_bot;
out.Ast2_top  = Ast2_top;
out.Ast2_bot  = Ast2_bot;
out.Fc_top    = Fc_top;
out.Fc_bot    = Fc_bot;
out.Sc_top    = Sc_top;
out.Sc_bot    = Sc_bot;

%% Bye, bye!
return;

function [NDes1, NDes2, Fc] = calc_NDes_Fc(N11, N22, N12)
F11 = min(N11, N22);
F22 = max(N11, N22);
F12 = N12;

if F11 + abs(F12) >= 0
    NDes1_  = F11 + abs(F12);
    NDes2_  = F22 + abs(F12);
    Fc      = -2*abs(F12);                  % concrete compressive force
else
    NDes1_  = 0;
    NDes2_  = F22 + F12^2/abs(F11);
    Fc      = -abs(F11)*(1 + (F12/F11)^2);  % concrete compressive force
end

if N11 <= N22
    NDes1 = NDes1_;
    NDes2 = NDes2_;
else
    NDes1 = NDes2_;
    NDes2 = NDes1_;
end

% Reinforcement forces must not be negative (force per unit length)
NDes1 = max(NDes1, 0);
NDes2 = max(NDes2, 0);

return
