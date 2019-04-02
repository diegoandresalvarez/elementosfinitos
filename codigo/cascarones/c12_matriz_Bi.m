%% Calculo de la matriz Bi

%% METODO 1
clear
clc
syms ti zi ui vi wi t1i t2i Ni dNi_dxi dNi_deta
invJ = sym('invJ_%d%d',[3 3]);
Ci   = sym('Ci_%d%d',[3 2]);

duvw_dxietazeta = [dNi_dxi *[eye(3)   -zi*Ci]  *[ui; vi; wi; t1i; t2i] ...
                   dNi_deta*[eye(3)   -zi*Ci]  *[ui; vi; wi; t1i; t2i] ...
                         Ni*[zeros(3) -ti*Ci/2]*[ui; vi; wi; t1i; t2i]].';
                      
duvw_dxyz = invJ*duvw_dxietazeta;

du_dx = duvw_dxyz(1,1);  dv_dx = duvw_dxyz(1,2);  dw_dx = duvw_dxyz(1,3);
du_dy = duvw_dxyz(2,1);  dv_dy = duvw_dxyz(2,2);  dw_dy = duvw_dxyz(2,3);
du_dz = duvw_dxyz(3,1);  dv_dz = duvw_dxyz(3,2);  dw_dz = duvw_dxyz(3,3);

ex = du_dx;    gxy = du_dy + dv_dx;
ey = dv_dy;    gxz = du_dz + dw_dx;
ez = dw_dz;    gyz = dv_dz + dw_dy;

def = expand([ ex; ey; ez; gxy; gxz; gyz ]);

ee = mat2cell(eye(5),[ 1 1 1 1 1 ], 5); % produce {[1 0 0 0 0];
                                        %          [0 1 0 0 0];
                                        %          ... 
                                        %          [0 0 0 0 1]}

Bi = subs(def, {ui; vi; wi; t1i; t2i}, ee)


%% METODO 2 (Comprobacion ecuacion caja 10.1 Onate)
syms invJ_11 invJ_21 invJ_31 invJ_12 invJ_22 invJ_32 invJ_13 invJ_23 invJ_33
Ni_1 = dNi_dxi*invJ_11 + dNi_deta*invJ_12;
Ni_2 = dNi_dxi*invJ_21 + dNi_deta*invJ_22;
Ni_3 = dNi_dxi*invJ_31 + dNi_deta*invJ_32;
Gi_1 = -(ti*invJ_13*Ni/2 + zi*Ni_1)*Ci;
Gi_2 = -(ti*invJ_23*Ni/2 + zi*Ni_2)*Ci;
Gi_3 = -(ti*invJ_33*Ni/2 + zi*Ni_3)*Ci;

Bi2 = ...
[ Ni_1,    0,    0,               Gi_1(1,1),               Gi_1(1,2)
     0, Ni_2,    0,               Gi_2(2,1),               Gi_2(2,2)
     0,    0, Ni_3,               Gi_3(3,1),               Gi_3(3,2)
  Ni_2, Ni_1,    0,   Gi_2(1,1) + Gi_1(2,1),   Gi_2(1,2) + Gi_1(2,2)
  Ni_3,    0, Ni_1,   Gi_3(1,1) + Gi_1(3,1),   Gi_3(1,2) + Gi_1(3,2)     
     0, Ni_3, Ni_2,   Gi_3(2,1) + Gi_2(3,1),   Gi_3(2,2) + Gi_2(3,2) ];

disp('Si da una matriz de ceros es porque Bi2 y Bi son iguales')  
disp(simple(Bi - Bi2))

% CONCLUSION: 
% ¡La matriz que aparece en el libro de Oñate está bien calculada y no tiene 
% errores! ¡esto es una gran sorpresa! :-P 
