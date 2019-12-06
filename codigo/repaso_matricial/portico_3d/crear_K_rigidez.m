clear, clc

syms AE EIy EIz GJ Iy Iz L L2 L3 
   
Kax = [ AE/L   -AE/L
       -AE/L    AE/L ];
 
Kv_tz = [ 12*EIz/L3   6*EIz/L2   -12*EIz/L3   6*EIz/L2
           6*EIz/L2   4*EIz/L     -6*EIz/L2   2*EIz/L
         -12*EIz/L3  -6*EIz/L2    12*EIz/L3  -6*EIz/L2
           6*EIz/L2   2*EIz/L     -6*EIz/L2   4*EIz/L   ];

Kw_ty = [ 12*EIy/L3   6*EIy/L2   -12*EIy/L3   6*EIy/L2
           6*EIy/L2   4*EIy/L     -6*EIy/L2   2*EIy/L
         -12*EIy/L3  -6*EIy/L2    12*EIy/L3  -6*EIy/L2
           6*EIy/L2   2*EIy/L     -6*EIy/L2   4*EIy/L   ];

Ktx = [ GJ/L   -GJ/L
       -GJ/L    GJ/L ];

K = sym(zeros(12,12));

idx = [1 7];
K(idx,idx) = K(idx,idx) + Kax;

idx = [2 6 8 12];
K(idx,idx) = K(idx,idx) + Kv_tz;

idx = [3 5 9 11];
K(idx,idx) = K(idx,idx) + Kw_ty;

idx = [4 10];
K(idx,idx) = K(idx,idx) + Ktx
