%% Calculo de la transformacion de deformaciones a otro sistema de coordenadas
syms exx eyy ezz gxy gxz gyz
exy = gxy/2;   exz = gxz/2;   eyz = gyz/2;

% syms a1 a2 a3 b1 b2 b3 g1 g2 g3     % como en el main.pdf
% T = [ a1 a2 a3
%       b1 b2 b3
%       g1 g2 g3 ]

syms lx ly lz mx my mz nx ny nz       % como en el libro de Onate
T = [ lx mx nx
      ly my ny
      lz mz nz ];
   
e = [ exx exy exz
      exy eyy eyz
      exz eyz ezz ];
      
ep = T.'*e*T;

expxp = ep(1,1);   expyp = ep(1,2);   expzp = ep(1,3);
                   eypyp = ep(2,2);   eypzp = ep(2,3);
                                      ezpzp = ep(3,3);

ep = expand([   expxp
                eypyp
                ezpzp
              2*expyp
              2*expzp
              2*eypzp ]);
                                      
%% Se calcula la matriz Q
ee = mat2cell(eye(6),[ 1 1 1 1 1 1 ], 6); % produce {[1 0 0 0 0 0];
                                          %          [0 1 0 0 0 0];
                                          %          ... 
                                          %          [0 0 0 0 0 1]}
                                        
Q = subs(ep, {exx; eyy; ezz; gxy; gxz; gyz}, ee)

