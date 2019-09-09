function [Bsp, det_J] = Bsp_lam_deg(e, xi, eta, zeta)
%% Calcula la matriz Bsp para el elemento de lamina degenerada en el punto 
%% de coordenadas (xi,eta,zeta)
%
% [Bsp, det_J] = Bsp_lam_deg(e, xi, eta, zeta)
% 
% e                      numero del EF
% (xi, eta, zeta)        punto de integracion de GL

%% Se calcula la matriz B
[B, Je, det_J] = B_lam_deg(e, xi, eta, zeta);

%% Se calcula x' y y'
xp = Je(1,:); % = dr_dxi  = [ dx_dxi;  dy_dxi;  dz_dxi  ]
yp = Je(2,:); % = dr_deta = [ dx_deta; dy_deta; dz_deta ]

%% Se calculan los elementos de la matriz T: lg, mn, ng
a1 = xp/norm(xp);
a2 = yp/norm(yp);

lg = a1;
ng = cross(a1, a2); ng = ng/norm(ng);
mg = cross(ng, lg);

%% Se calcula la matriz Q2
lx = lg(1);   mx = mg(1);   nx = ng(1);
ly = lg(2);   my = mg(2);   ny = ng(2);
lz = lg(3);   mz = mg(3);   nz = ng(3);

% la matriz Q se de dedujo con el programa c12_matriz_Q.m (eliminando la
% tercera fila)
Q2 = [ ... 
%   lx^2,    ly^2,    lz^2,         lx*ly,         lx*lz,         ly*lz
%   mx^2,    my^2,    mz^2,         mx*my,         mx*mz,         my*mz
%2*lx*mx, 2*ly*my, 2*lz*mz, lx*my + ly*mx, lx*mz + lz*mx, ly*mz + lz*my ];
 2*lx*nx, 2*ly*ny, 2*lz*nz, lx*ny + ly*nx, lx*nz + lz*nx, ly*nz + lz*ny
 2*mx*nx, 2*my*ny, 2*mz*nz, mx*ny + my*nx, mx*nz + mz*nx, my*nz + mz*ny ];

%% Se calcula la matriz Bpp
Bsp = Q2*B;

return;
