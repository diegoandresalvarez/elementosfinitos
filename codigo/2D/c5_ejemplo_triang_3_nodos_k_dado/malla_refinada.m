% Malla elaborada por David Felipe Cano Perdomo

xnod = [ 0.000     0.000 % posicion de los nodos
         0.050     0.000    
         0.100     0.000 % la fila representa el numero del nodo
         0.150     0.000 % la columna representa la coordenada X=1 o Y=2
         0.200     0.000
         0.250     0.000
         0.300     0.000
         0.350     0.000
         0.400     0.000
         0.450     0.000
         0.500     0.000
         0.550     0.000
         0.600     0.000
         0.650     0.000
         0.700     0.000
         0.750     0.000
         0.800     0.000
         0.000     0.035
         0.050     0.035
         0.100     0.035
         0.150     0.035
         0.200     0.035
         0.250     0.035
         0.300     0.035
         0.350     0.035
         0.400     0.035
         0.450     0.035
         0.500     0.035
         0.550     0.035
         0.600     0.035
         0.650     0.035
         0.700     0.035
         0.750     0.035
         0.800     0.035
         0.000     0.070
         0.050     0.070
         0.100     0.070
         0.150     0.070
         0.200     0.070
         0.250     0.070
         0.300     0.070
         0.350     0.070
         0.400     0.070
         0.450     0.070
         0.500     0.070
         0.550     0.070
         0.600     0.070
         0.650     0.070
         0.700     0.070
         0.750     0.070
         0.800     0.070
         0.000     0.105
         0.050     0.105
         0.100     0.105
         0.150     0.105
         0.200     0.105
         0.250     0.105
         0.300     0.105
         0.350     0.105
         0.400     0.105
         0.450     0.105
         0.500     0.105
         0.550     0.105
         0.600     0.105
         0.650     0.105
         0.700     0.105
         0.750     0.105
         0.800     0.105
         0.000     0.140
         0.050     0.140
         0.100     0.140
         0.150     0.140
         0.200     0.140
         0.250     0.140
         0.300     0.140
         0.350     0.140
         0.400     0.140
         0.450     0.140
         0.500     0.140
         0.550     0.140
         0.600     0.140
         0.650     0.140
         0.700     0.140
         0.750     0.140
         0.800     0.140
         0.000     0.175
         0.050     0.175
         0.100     0.175
         0.150     0.175
         0.200     0.175
         0.250     0.175
         0.300     0.175
         0.350     0.175
         0.400     0.175
         0.450     0.175
         0.500     0.175
         0.550     0.175
         0.600     0.175
         0.650     0.175
         0.700     0.175
         0.750     0.175
         0.800     0.175
         0.000     0.210
         0.050     0.210
         0.100     0.210
         0.150     0.210
         0.200     0.210
         0.250     0.210
         0.300     0.210
         0.350     0.210
         0.400     0.210
         0.450     0.210
         0.500     0.210
         0.550     0.210
         0.600     0.210
         0.650     0.210
         0.700     0.210
         0.750     0.210
         0.800     0.210
         0.000     0.245
         0.050     0.245
         0.100     0.245
         0.150     0.245
         0.200     0.245
         0.250     0.245
         0.300     0.245
         0.350     0.245
         0.400     0.245
         0.450     0.245
         0.500     0.245
         0.550     0.245
         0.600     0.245
         0.650     0.245
         0.700     0.245
         0.750     0.245
         0.800     0.245
         0.000     0.280
         0.050     0.280
         0.100     0.280
         0.150     0.280
         0.200     0.280
         0.250     0.280
         0.300     0.280
         0.350     0.280
         0.400     0.280
         0.450     0.280
         0.500     0.280
         0.550     0.280
         0.600     0.280
         0.650     0.280
         0.700     0.280
         0.750     0.280
         0.800     0.280 ];

nno = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 2*nno;       % numero de grados de libertad (dos por nodo)
gdl = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad

LaG = [   1   2  19 % definicion de elementos finitos con respecto a nodos
         19  18   1
         18  19  36 % la fila representa el numero del elemento
         36  35  18 % la columna representa el numero de nodo local 1,2 o 3
         35  36  53
         53  52  35
         52  53  70
         70  69  52
         69  70  87
         87  86  69
         86  87 104
        104 103  86
        103 104 121
        121 120 103
        120 121 138
        138 137 120
          2   3  20
         20  19   2
         19  20  37
         37  36  19
         36  37  54
         54  53  36
         53  54  71
         71  70  53
         70  71  88
         88  87  70
         87  88 105
        105 104  87
        104 105 122
        122 121 104
        121 122 139
        139 138 121
          3   4  21
         21  20   3
         20  21  38
         38  37  20
         37  38  55
         55  54  37
         54  55  72
         72  71  54
         71  72  89
         89  88  71
         88  89 106
        106 105  88
        105 106 123
        123 122 105
        122 123 140
        140 139 122
          4   5  22
         22  21   4
         21  22  39
         39  38  21
         38  39  56
         56  55  38
         55  56  73
         73  72  55
         72  73  90
         90  89  72
         89  90 107
        107 106  89
        106 107 124
        124 123 106
        123 124 141
        141 140 123
          5   6  23
         23  22   5
         22  23  40
         40  39  22
         39  40  57
         57  56  39
         56  57  74
         74  73  56
         73  74  91
         91  90  73
         90  91 108
        108 107  90
        107 108 125
        125 124 107
        124 125 142
        142 141 124
          6   7  24
         24  23   6
         23  24  41
         41  40  23
         40  41  58
         58  57  40
         57  58  75
         75  74  57
         74  75  92
         92  91  74
         91  92 109
        109 108  91
        108 109 126
        126 125 108
        125 126 143
        143 142 125
          7   8  25
         25  24   7
         24  25  42
         42  41  24
         41  42  59
         59  58  41
         58  59  76
         76  75  58
         75  76  93
         93  92  75
         92  93 110
        110 109  92
        109 110 127
        127 126 109
        126 127 144
        144 143 126
          8   9  26
         26  25   8
         25  26  43
         43  42  25
         42  43  60
         60  59  42
         59  60  77
         77  76  59
         76  77  94
         94  93  76
         93  94 111
        111 110  93
        110 111 128
        128 127 110
        127 128 145
        145 144 127
          9  10  27
         27  26   9
         26  27  44
         44  43  26
         43  44  61
         61  60  43
         60  61  78
         78  77  60
         77  78  95
         95  94  77
         94  95 112
        112 111  94
        111 112 129
        129 128 111
        128 129 146
        146 145 128
         10  11  28
         28  27  10
         27  28  45
         45  44  27
         44  45  62
         62  61  44
         61  62  79
         79  78  61
         78  79  96
         96  95  78
         95  96 113
        113 112  95
        112 113 130
        130 129 112
        129 130 147
        147 146 129
         11  12  29
         29  28  11
         28  29  46
         46  45  28
         45  46  63
         63  62  45
         62  63  80
         80  79  62
         79  80  97
         97  96  79
         96  97 114
        114 113  96
        113 114 131
        131 130 113
        130 131 148
        148 147 130
         12  13  30
         30  29  12
         29  30  47
         47  46  29
         46  47  64
         64  63  46
         63  64  81
         81  80  63
         80  81  98
         98  97  80
         97  98 115
        115 114  97
        114 115 132
        132 131 114
        131 132 149
        149 148 131
         13  14  31
         31  30  13
         30  31  48
         48  47  30
         47  48  65
         65  64  47
         64  65  82
         82  81  64
         81  82  99
         99  98  81
         98  99 116
        116 115  98
        115 116 133
        133 132 115
        132 133 150
        150 149 132
         14  15  32
         32  31  14
         31  32  49
         49  48  31
         48  49  66
         66  65  48
         65  66  83
         83  82  65
         82  83 100
        100  99  82
         99 100 117
        117 116  99
        116 117 134
        134 133 116
        133 134 151
        151 150 133
         15  16  33
         33  32  15
         32  33  50
         50  49  32
         49  50  67
         67  66  49
         66  67  84
         84  83  66
         83  84 101
        101 100  83
        100 101 118
        118 117 100
        117 118 135
        135 134 117
        134 135 152
        152 151 134
         16  17  34
         34  33  16
         33  34  51
         51  50  33
         50  51  68
         68  67  50
         67  68  85
         85  84  67
         84  85 102
        102 101  84
        101 102 119
        119 118 101
        118 119 136
        136 135 118
        135 136 153
        153 152 135];

nef = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se definen las restricciones 
%             gdl       desplazamiento(m)
restric = [   gdl(1,X)          0
              gdl(1,Y)          0
              gdl(17,Y)          0  ];
           
%% Se definen las cargas distribuidas 
%             [ elemento  lado  tix   tiy  tjx  tjy  ]
carga_distr = [     16     12   0   -9.375  0   0
                    32     12   0   -8.75   0   -9.375
                    48     12   0   -8.125  0   -8.75
                    64     12   0   -7.5    0   -8.125
                    80     12   0   -6.875  0   -7.5
                    96     12   0   -6.25   0   -6.875
                    112    12   0   -5.625  0   -6.25
                    128    12   0   -5      0   -5.625];
                                 
nlcd = size(carga_distr,1); % numero de lados con carga distribuida
         
%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
f(gdl(149,X)) = 7000*cosd(60);  % carga puntual en el nodo 14 dir X
f(gdl(149,Y)) = 7000*sind(60);  % carga puntual en el nodo 14 dir Y        
    
