syms xi eta

% Funciones de forma del EF de 4 nodos
N = sym(zeros(1,4));
N(1) = (1 - xi)*(1 - eta)/4;
N(2) = (1 + xi)*(1 - eta)/4;
N(3) = (1 + xi)*(1 + eta)/4;
N(4) = (1 - xi)*(1 + eta)/4;

% Se definen los vectores de desplazamientos y giros nodales
we  = sym('w',  [4 1]);   w  = N*we;
txe = sym('tx', [4 1]);   tx = N*txe;
tye = sym('ty', [4 1]);   ty = N*tye;

% Se calculan las deformaciones gxz y se extraen sus coeficientes
L = 2;
dxi_dx = 2/L;    deta_dx = 0;
dxi_dy = 0;      deta_dy = 2/L;
dw_dx = diff(w,xi)*dxi_dx + diff(w,eta)*deta_dx;
dw_dy = diff(w,xi)*dxi_dy + diff(w,eta)*deta_dy;
gxz = dw_dx - tx;
gyz = dw_dy - ty;

% Se evaluan las constantes alpha y alphab
alpha = sym(zeros(4,1));
alpha(1) = feval(symengine, 'coeff', gxz, '[xi,eta]', '[0,0]');
alpha(2) = feval(symengine, 'coeff', gxz, '[xi,eta]', '[1,0]');
alpha(3) = feval(symengine, 'coeff', gxz, '[xi,eta]', '[0,1]');
alpha(4) = feval(symengine, 'coeff', gxz, '[xi,eta]', '[1,1]');

alphab = sym(zeros(4,1));
alphab(1) = feval(symengine, 'coeff', gyz, '[xi,eta]', '[0,0]');
alphab(2) = feval(symengine, 'coeff', gyz, '[xi,eta]', '[1,0]');
alphab(3) = feval(symengine, 'coeff', gyz, '[xi,eta]', '[0,1]');
alphab(4) = feval(symengine, 'coeff', gyz, '[xi,eta]', '[1,1]');
% expand(gxz - (alpha(1) + alpha(2)*xi + alpha(3)*eta + alpha(4)*xi*eta))
% expand(gyz - (alphab(1) + alphab(2)*xi + alphab(3)*eta + alphab(4)*xi*eta))

cond_xz = simplify(alpha  == 0)
cond_yz = simplify(alphab == 0)
cond_xz(2) - cond_xz(4)
cond_yz(3) - cond_yz(4)

