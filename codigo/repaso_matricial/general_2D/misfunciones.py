# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# %% se definen algunas constantes
X,  Y = 0,  1
X1, Y1, M1, X2, Y2, M2 = range(6)

# %%
def calc_Te(x1,y1, x2,y2):
    '''
    Matriz de transformación T para un EF de pórtico/cercha con 3 gdl por nodo

    PARAMETROS:
    (x1,y1) y (x2,y2) son las coordenadas de los extremos de la barra
    '''
    L = np.hypot(x2-x1, y2-y1)    # longitud de la barra
    c = (x2-x1)/L; s = (y2-y1)/L  # coseno y seno de la inclinación
    Te = np.array([[ c,  s,  0,  0,  0,  0],
                   [-s,  c,  0,  0,  0,  0],
                   [ 0,  0,  1,  0,  0,  0],
                   [ 0,  0,  0,  c,  s,  0],
                   [ 0,  0,  0, -s,  c,  0],
                   [ 0,  0,  0,  0,  0,  1]])

    return Te

# %%
def calc_feloc(tipo, L, b1,b2, q1,q2):
    '''
    Calcula el vector de fuerzas nodales equivalentes de un EF de pórtico o de 
    cercha en coordenadas locales

    PARAMETROS:
    tipo    = 'EE' (pórtico)
            = 'RR' (cercha)
            = 'RE' (rotula/empotrado)
            = 'ER' (empotrado/rotula)
    L       longitud de la barra
    b1, b2  carga axial    en x1 y x2 (especificada en coordenadas locales)
    q1, q2  carga vertical en x1 y x2 (especificada en coordenadas locales)
    '''
    # %% se calculan las fuerzas nodales de equilibrio en coordenadas locales
    if tipo in ['EE', 'ER_PH', 'ER_HP']:
        X1 =  (L*(2*b1 + b2))/6
        Y1 =  (L*(7*q1 + 3*q2))/20
        M1 =  (L**2 * (3*q1 + 2*q2))/60
        X2 =  (L*(b1 + 2*b2))/6
        Y2 =  (L*(3*q1 + 7*q2))/20
        M2 = -(L**2 * (2*q1 + 3*q2))/60
    elif tipo == 'RR':
        X1 = (L*(2*b1 + b2))/6
        Y1 = (L*(2*q1 + q2))/6
        M1 = 0
        X2 = (L*(b1 + 2*b2))/6 
        Y2 = (L*(q1 + 2*q2))/6
        M2 = 0        
    elif tipo == 'RE':
        X1 = (L*(2*b1 + b2))/6
        Y1 = (L*(11*q1 + 4*q2))/40
        M1 = 0
        X2 = (L*(b1 + 2*b2))/6
        Y2 = (L*(9*q1 + 16*q2))/40
        M2 = -(L**2*(7*q1 + 8*q2))/120        
    elif tipo == 'ER':
        X1 = (L*(2*b1 + b2))/6
        Y1 = (L*(16*q1 + 9*q2))/40
        M1 = (L**2*(8*q1 + 7*q2))/120
        X2 = (L*(b1 + 2*b2))/6
        Y2 = (L*(4*q1 + 11*q2))/40
        M2 = 0        
    else:       
        raise ValueError("Tipo de EF no soportado")
    
    feloc = np.array([ X1, Y1, M1, X2, Y2, M2 ]) 
    
    return feloc

# %%
def calc_Keloc(tipo, L, A, E, I):
    '''
    Calcula la matriz de rigidez local de un EF de pórtico o de cercha en 
    coordenadas locales

    PARAMETROS:
    tipo    = 'EE' (pórtico)
            = 'RR' (cercha)
            = 'RE' (rotula/empotrado)
            = 'ER' (empotrado/rotula)
    A      área  
    E      módulo de elasticidad
    I      inercia
    L      longitud de la barra
    '''       
    
    AE = A*E;       L2=L**2
    EI = E*I;       L3=L**3
    
    # matriz de rigidez local expresada en el sistema de coordenadas locales
    if tipo in ['EE', 'ER_PH', 'ER_HP']:
        Keloc = np.array([
             [ AE/L,   0      ,   0      ,  -AE/L,    0      ,   0      ],  
             [ 0   ,  12*EI/L3,   6*EI/L2,   0   ,  -12*EI/L3,   6*EI/L2],
             [ 0   ,   6*EI/L2,   4*EI/L ,   0   ,   -6*EI/L2,   2*EI/L ],
             [-AE/L,   0      ,   0      ,   AE/L,    0      ,   0      ],
             [ 0   , -12*EI/L3,  -6*EI/L2,   0   ,   12*EI/L3,  -6*EI/L2],
             [ 0   ,   6*EI/L2,   2*EI/L ,   0   ,   -6*EI/L2,   4*EI/L ]])
    elif tipo == 'RR':
        k = AE/L
        Keloc = np.array([[ k,  0,  0, -k,  0,  0],
                          [ 0,  0,  0,  0,  0,  0],
                          [ 0,  0,  0,  0,  0,  0],                          
                          [-k,  0,  0,  k,  0,  0],
                          [ 0,  0,  0,  0,  0,  0],                          
                          [ 0,  0,  0,  0,  0,  0]])        
    elif tipo == 'ER':
        Keloc = np.array([
             [ AE/L,   0      ,   0      ,  -AE/L,    0      ,   0      ],  
             [ 0   ,   3*EI/L3,   3*EI/L2,   0   ,   -3*EI/L3,   0      ],
             [ 0   ,   3*EI/L2,   3*EI/L ,   0   ,   -3*EI/L2,   0      ],
             [-AE/L,   0      ,   0      ,   AE/L,    0      ,   0      ],
             [ 0   ,  -3*EI/L3,  -3*EI/L2,   0   ,    3*EI/L3,   0      ],
             [ 0   ,   0      ,   0      ,   0   ,    0      ,   0      ]])
    elif tipo == 'RE':
        Keloc = np.array([
             [ AE/L,   0      ,   0      ,  -AE/L,    0      ,   0      ],  
             [ 0   ,   3*EI/L3,   0      ,   0   ,   -3*EI/L3,   3*EI/L2],
             [ 0   ,   0      ,   0      ,   0   ,    0      ,   0      ],
             [-AE/L,   0      ,   0      ,   AE/L,    0      ,   0      ],
             [ 0   ,  -3*EI/L3,   0      ,   0   ,    3*EI/L3,  -3*EI/L2],
             [ 0   ,   3*EI/L2,   0      ,   0   ,   -3*EI/L2,   3*EI/L ]])        
    else:
        raise ValueError("Tipo de EF no soportado")

    return Keloc

#%%
def crear_C_g_enlace_rigido(xnod, idx_ER):
    '''
    Esta función retorna las matrices C y g asociadas al enlace rígido.

    Para tal fin implementa dentro de C y g las ecuaciones:

    u(h) = u(p) - t(p)*dy
    v(h) = v(p) + t(p)*dx
    t(h) = t(p)

    USO:
    C, g = crear_C_g_enlace_rigido(xnod, idx_ER)    

    PARAMETROS DE ENTRADA:
    xnod   = matriz de nnod x 2 con las coordenadas x,y de cada nodo
    idx_ER = matriz de num_ER x 2: columna 1=nodo padre, columna 2=nodo hijo

    PARAMETROS DE SALIDA:
    C
    g
    '''

    # extrae el número de enlaces rígidos:
    num_ER = idx_ER.shape[0]

    # si no existen enlaces rígidos, sálgase
    if num_ER == 0: return None, None, None
    
    # 
    nno  = xnod.shape[0] # número de nodos
    ngdl = 3*nno         # número de grados de libertad por nodo = [X, Y, TH]
   
    # matrices que sirven para definir las restricciones
    C = np.zeros((3*num_ER, ngdl))
    g = np.zeros((3*num_ER))
    
    gdl_hijos = []
    
    for i in range(num_ER):
        nodo_padre, nodo_hijo  = idx_ER[i,:]
        
        dx = xnod[nodo_hijo, X] - xnod[nodo_padre, X]
        dy = xnod[nodo_hijo, Y] - xnod[nodo_padre, Y]
        
        # se calculan los idx de los GDL del nodo padre y del nodo hijo
        gdl_up = 3*nodo_padre;           gdl_uh = 3*nodo_hijo    
        gdl_vp = 3*nodo_padre + 1;       gdl_vh = 3*nodo_hijo + 1
        gdl_tp = 3*nodo_padre + 2;       gdl_th = 3*nodo_hijo + 2
        gdl_hijos.extend([gdl_uh, gdl_vh, gdl_th])
        
        # se escriben las ecuaciones asociadas a los enlaces rígidos
        C[3*i+0, [gdl_uh, gdl_up, gdl_tp]] = [1, -1, +dy] # u(h) - u(p) + t(p)*dy == 0
        C[3*i+1, [gdl_vh, gdl_vp, gdl_tp]] = [1, -1, -dx] # v(h) - v(p) - t(p)*dx == 0
        C[3*i+2, [gdl_th, gdl_tp]]         = [1, -1]      # t(h) - t(p)           == 0
    
    return C, g, gdl_hijos

# %%
def resolver_Ka_f_q(K, f, c, ac, C=None, g=None, gdl_hijos=None):
    '''
    Esta función resuelve el sistema de ecuaciones Ka - f = q

    PARAMETROS DE ENTRADA:
    K  = matriz de rigidez global
    f  = vector de fuerzas nodales equivalentes    
    c  = GDL del desplazamiento conocidos (GDL de los apoyos)
    ac = desplazamientos conocidos
    C  = 
    g  =
    gdl_hijos = índices de los GDL de los hijos


    PARAMETROS DE SALIDA:
    a  = desplazamientos
    q  = vector de fuerzas nodales de equilibrio del elemento
    '''

    # número de grados de libertad
    ngdl = f.shape[0]

    #| qd |   | Kcc Kcd || ac |   | fd |    Recuerde que siempre qc=0
    #|    | = |         ||    | - |    |
    #| qc |   | Kdc Kdd || ad |   | fc |

    if C is None: 
        #%% si no hay restricciones/enlaces rígidos:

        # grados de libertad del desplazamiento desconocidos
        d = np.setdiff1d(np.arange(ngdl), c)

        # extraigo las submatrices y especifico las cantidades conocidas
        Kcc = K[np.ix_(c,c)];      Kcd = K[np.ix_(c,d)];      fd = f[c]
        Kdc = K[np.ix_(d,c)];      Kdd = K[np.ix_(d,d)];      fc = f[d]
    
        # resuelvo el sistema de ecuaciones
        ad = np.linalg.solve(Kdd, fc - Kdc@ac) # desplazamientos desconocidos
        qd = Kcc@ac + Kcd@ad - fd              # fuerzas de equilibrio desconocidas
    
        # armo los vectores de desplazamientos (a) y fuerzas (q)
        a = np.zeros(ngdl);   a[c] = ac;   a[d] = ad # desplazamientos
        q = np.zeros(ngdl);   q[c] = qd              # fuerzas nodales equivalentes
    else:
        # %%si hay restricciones/enlaces rígidos, use el método 1: 
        # transformacion de K*a-f = q (ver diapositivas 19)

        e = gdl_hijos # idx de los GDL asociados a los nodos hijos
        
        # se calculan los idx de los GDL asociados a los nodos a retener
        # (incluye los nodos maestros)
        r = np.setdiff1d(np.arange(ngdl), e).tolist()
        
        # se verifica que ningún GDL de desplazamiento conocido sea un GDL hijo
        if np.any(np.intersect1d(c, e)):
            raise Exception('No pueden haber GDL del empotramiento en los GDL hijos')
        
        num_hijos = C.shape[0]         # número de GDL hijos     (c)
        n_c       = ngdl - num_hijos   # número de GDL a retener (n-c)
       
        # se extraen las submatrices Cr y Ce de la matriz de restricciones C
        Cr = C[:,r];   Ce = C[:,e]
        
        # Se verifica que Ce no sea singular
        if np.abs(np.linalg.det(Ce)) < 1e-5:
            raise Exception('Ce debe ser invertible')
        
        I  = np.eye(n_c)
        T  = np.r_[  I,
                    -np.linalg.solve(Ce, Cr) ] # = -inv(Ce)*Cr
        
        O  = np.zeros(n_c)
        g0 = np.r_[  O,
                     np.linalg.solve(Ce, g)  ] # =  inv(Ce)*g
        
        # se reordenan los GDL de K y de f
        idx_re = np.r_[r, e].tolist()
        
        Kre = K[np.ix_(idx_re, idx_re)]
        fre = f[idx_re]
           
        Kr = T.T @ Kre @ T
        fr = T.T @ (fre - Kre@g0)
        
        # Se definen GDL conocidos y desconocidos asociados a los desplazamientos
        # idx con la numeración original (del grafico)
        # observe que en Kr*ar - fr = qr no hay nodos esclavos
        d = np.setdiff1d(r, c)
        
        # nuevos idx de c y d con la numeración del reordenamiento [r,e], esto es idx_re
        c_re = [ r.index(i) for i in c ]
        d_re = [ r.index(i) for i in d ]
        
        # Se descomponen los vectores a, f y la matriz K
        Krcc = Kr[np.ix_(c_re,c_re)];      Krcd = Kr[np.ix_(c_re,d_re)]
        Krdc = Kr[np.ix_(d_re,c_re)];      Krdd = Kr[np.ix_(d_re,d_re)]
        frc  = fr[d_re];                   frd  = fr[c_re]
        
        # se reordenan los GDL de a y se extrae arc
        a = np.zeros(ngdl)
        a[c] = ac
        ar   = a[r]
        arc  = ar[c_re]
        
        # Se calculan los vectores ard y qrd (de los GDL a retener)
        # resuelvo el sistema de ecuaciones
        ard = np.linalg.solve(Krdd, frc - Krdc@arc) # desplazamientos desconocidos
        qrd = Krcc@arc + Krcd@ard - frd             # fuerzas de equilibrio desconocidas
        
        # armo los vectores de desplazamientos (a_re) y fuerzas (q)
        ar = np.zeros(n_c)     # desplazamientos
        ar[c_re] = arc
        ar[d_re] = ard
        a_re     = T@ar + g0
        
        # se calcula el vector q_re
        q_re       = np.zeros(ngdl)  # fuerzas nodales equivalentes
        q_re[c_re] = qrd
        
        # se reordena a de la indexación re a la indexación original (del gráfico)
        idx_orig = [ idx_re.index(i) for i in range(ngdl) ]
        a = a_re[idx_orig]
        q = q_re[idx_orig]
        
    # %% bye, bye!!!
    return a, q

# %%
def dibujar_deformada(tipo, A, E, I, x1,y1, x2,y2, b1,b2, q1,q2,
                    qe, ae, esc_def, esc_faxial, esc_V, esc_M, con_lineas=True):
    '''
    Esta función dibuja el elemento de pórtico/cercha deformado junto con sus 
    respectivos diagramas de fuerza axial, fuerza cortante y momento flector.

    El diagrama de momento flector se grafica en el lado opuesto de la fibra
    a tracción

    PARAMETROS DE ENTRADA (junto con algunos ejemplos):
    tipo    = 'EE' (pórtico)
            = 'RR' (cercha)
            = 'RE' (rotula/empotrado)
            = 'ER' (empotrado/rotula)
    A = area
    E = E
    I = Ix local
    (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
    b1, b2    carga axial    en x1 y x2
    q1, q2    carga vertical en x1 y x2  
    qe = np.array(
         [ 0.01,        # U1, V1, M1 reacciones del nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,        # U2, V2, M2 reacciones del nodo 2 en coord. locales
           0.02,
          -0.07 ])
    ae = np.array(
         [ 0.01,        # u1, v1, t1 desplazamientos nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,        # u2, v2, t2 desplazamientos nodo 2 en coord. locales
           0.02,
          -0.07 ])
    esc_def    = 10     # escalamiento de la deformada
    esc_faxial = 10     # escalamiento del diagrama de axiales
    esc_V      = 10     # escalamiento del diagrama de cortantes
    esc_M      = 10     # escalamiento del diagrama de momentos
    con_lineas = {True/False} # graficar lineas verticales en diagramas de M/V/A
    '''
    npuntos = 1001
    L = np.hypot(x2-x1, y2-y1)
    x = s = np.linspace(0, L, npuntos)  
    
    # se asignan los desplazamientos en coordenadas locales
    u1, v1, t1, u2, v2, t2 = ae
    
    # se calculan mediante fórmulas las fuerzas, el momento y los desplazamientos
    #q = q1 - (x*(q1 - q2))/L
    if tipo in ['EE', 'ER_PH', 'ER_HP']:
        axial = ((b1 - b2)*x**2)/(2*L) - b1*x + (2*L**2*b1 + L**2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L)
        V     = q1*x - (7*L**4*q1 + 3*L**4*q2 - 240*E*I*v1 + 240*E*I*v2 - 120*E*I*L*t1 - 120*E*I*L*t2)/(20*L**3) - (x**2*(q1 - q2))/(2*L)
        M     = (3*L**5*q1 + 2*L**5*q2 - 10*L**2*q1*x**3 + 30*L**3*q1*x**2 + 10*L**2*q2*x**3 - 21*L**4*q1*x - 9*L**4*q2*x - 240*E*I*L**2*t1 - 120*E*I*L**2*t2 - 360*E*I*L*v1 + 360*E*I*L*v2 + 720*E*I*v1*x - 720*E*I*v2*x + 360*E*I*L*t1*x + 360*E*I*L*t2*x)/(60*L**3)
        u     = (b1*x**3 - b2*x**3 - 3*L*b1*x**2 + 2*L**2*b1*x + L**2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L)
        v     = (5*L**3*q1*x**4 - L**2*q1*x**5 - 7*L**4*q1*x**3 + 3*L**5*q1*x**2 + L**2*q2*x**5 - 3*L**4*q2*x**3 + 2*L**5*q2*x**2 + 120*E*I*L**3*v1 + 240*E*I*v1*x**3 - 240*E*I*v2*x**3 + 120*E*I*L*t1*x**3 + 120*E*I*L**3*t1*x + 120*E*I*L*t2*x**3 - 360*E*I*L*v1*x**2 + 360*E*I*L*v2*x**2 - 240*E*I*L**2*t1*x**2 - 120*E*I*L**2*t2*x**2)/(120*E*I*L**3)
        #t    = (20*L**3*q1*x**3 - 5*L**2*q1*x**4 - 21*L**4*q1*x**2 + 5*L**2*q2*x**4 - 9*L**4*q2*x**2 + 6*L**5*q1*x + 4*L**5*q2*x + 120*E*I*L**3*t1 + 720*E*I*v1*x**2 - 720*E*I*v2*x**2 - 720*E*I*L*v1*x + 720*E*I*L*v2*x + 360*E*I*L*t1*x**2 - 480*E*I*L**2*t1*x + 360*E*I*L*t2*x**2 - 240*E*I*L**2*t2*x)/(120*E*I*L**3)
    elif tipo == 'RR':
        axial = ((b1 - b2)*x**2)/(2*L) - b1*x + (2*L**2*b1 + L**2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L)
        V     = q1*x - (L*q2)/6 - (L*q1)/3 - (x**2*(q1 - q2))/(2*L)
        M     = -(x*(L - x)*(2*L*q1 + L*q2 - q1*x + q2*x))/(6*L)
        u     = (b1*x**3 - b2*x**3 - 3*L*b1*x**2 + 2*L**2*b1*x + L**2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L)
        v     = (3*q2*x**5 - 3*q1*x**5 - 20*L**2*q1*x**3 - 10*L**2*q2*x**3 + 15*L*q1*x**4 + 8*L**4*q1*x + 7*L**4*q2*x + 360*E*I*L*v1 - 360*E*I*v1*x + 360*E*I*v2*x)/(360*E*I*L)
        #t    = (8*L**4*q1 + 7*L**4*q2 - 15*q1*x**4 + 15*q2*x**4 - 60*L**2*q1*x**2 - 30*L**2*q2*x**2 - 360*E*I*v1 + 360*E*I*v2 + 60*L*q1*x**3)/(360*E*I*L)   
    elif tipo == 'RE':
        axial = ((b1 - b2)*x**2)/(2*L) - b1*x + (2*L**2*b1 + L**2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L)
        V     = q1*x - (11*L**4*q1 + 4*L**4*q2 - 120*E*I*v1 + 120*E*I*v2 - 120*E*I*L*t2)/(40*L**3) - (x**2*(q1 - q2))/(2*L)
        M     = -(x*(33*L**4*q1 + 12*L**4*q2 + 20*L**2*q1*x**2 - 20*L**2*q2*x**2 - 360*E*I*v1 + 360*E*I*v2 - 60*L**3*q1*x - 360*E*I*L*t2))/(120*L**3)
        u     = (b1*x**3 - b2*x**3 - 3*L*b1*x**2 + 2*L**2*b1*x + L**2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L)
        v     = (10*L**3*q1*x**4 - 2*L**2*q1*x**5 - 11*L**4*q1*x**3 + 2*L**2*q2*x**5 - 4*L**4*q2*x**3 + 3*L**6*q1*x + 2*L**6*q2*x + 240*E*I*L**3*v1 + 120*E*I*v1*x**3 - 120*E*I*v2*x**3 + 120*E*I*L*t2*x**3 - 120*E*I*L**3*t2*x - 360*E*I*L**2*v1*x + 360*E*I*L**2*v2*x)/(240*E*I*L**3)
        #t    = (3*L**6*q1 + 2*L**6*q2 - 10*L**2*q1*x**4 + 40*L**3*q1*x**3 - 33*L**4*q1*x**2 + 10*L**2*q2*x**4 - 12*L**4*q2*x**2 - 120*E*I*L**3*t2 - 360*E*I*L**2*v1 + 360*E*I*L**2*v2 + 360*E*I*v1*x**2 - 360*E*I*v2*x**2 + 360*E*I*L*t2*x**2)/(240*E*I*L**3)
    elif tipo == 'ER':
        axial = ((b1 - b2)*x**2)/(2*L) - b1*x + (2*L**2*b1 + L**2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L)
        V     = q1*x - (16*L**4*q1 + 9*L**4*q2 - 120*E*I*v1 + 120*E*I*v2 - 120*E*I*L*t1)/(40*L**3) - (x**2*(q1 - q2))/(2*L)
        M     = -((L - x)*(20*L**2*q2*x**2 - 7*L**4*q2 - 20*L**2*q1*x**2 - 8*L**4*q1 + 360*E*I*v1 - 360*E*I*v2 + 40*L**3*q1*x + 20*L**3*q2*x + 360*E*I*L*t1))/(120*L**3)
        u     = (b1*x**3 - b2*x**3 - 3*L*b1*x**2 + 2*L**2*b1*x + L**2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L)
        v     = (10*L**3*q1*x**4 - 2*L**2*q1*x**5 - 16*L**4*q1*x**3 + 8*L**5*q1*x**2 + 2*L**2*q2*x**5 - 9*L**4*q2*x**3 + 7*L**5*q2*x**2 + 240*E*I*L**3*v1 + 120*E*I*v1*x**3 - 120*E*I*v2*x**3 + 120*E*I*L*t1*x**3 + 240*E*I*L**3*t1*x - 360*E*I*L*v1*x**2 + 360*E*I*L*v2*x**2 - 360*E*I*L**2*t1*x**2)/(240*E*I*L**3)
        #t    = (40*L**3*q1*x**3 - 10*L**2*q1*x**4 - 48*L**4*q1*x**2 + 10*L**2*q2*x**4 - 27*L**4*q2*x**2 + 16*L**5*q1*x + 14*L**5*q2*x + 240*E*I*L**3*t1 + 360*E*I*v1*x**2 - 360*E*I*v2*x**2 - 720*E*I*L*v1*x + 720*E*I*L*v2*x + 360*E*I*L*t1*x**2 - 720*E*I*L**2*t1*x)/(240*E*I*L**3)
    else:
        raise Exception('Tipo de EF no soportado')
    #theta = np.arctan(t)  # Angulo de giro [rad]
    
    # rotación de la solución antes de dibujar
    cos_ang = (x2-x1)/L
    sin_ang = (y2-y1)/L
    T       = np.array([[cos_ang,  -sin_ang],  # matriz de rotación
                        [sin_ang,   cos_ang]]) 
      
    #%% Dibujar de deformada
    plt.figure(2)
    pos = T @ np.array([s + esc_def*u,
                            esc_def*v])
    
    ss = pos[X,:] + x1
    yy = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-', ss, yy, 'r-', linewidth=2)
    
    #%% Dibujar los diagramas de fuerza axial 
    plt.figure(3)
    pos = T @ np.array([s               ,
                        esc_faxial*axial]) # escalamiento del diagrama
    
    ss = pos[X,:] + x1
    aa = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-')
    if con_lineas:
        plt.plot(np.r_[x1, ss, x2], np.r_[y1, aa, y2], 'r-', linewidth=2)
    else:
        plt.plot(ss, aa, 'r-', linewidth=2)

    plt.text(ss[ 0], aa[ 0], f'{-qe[X1]:.4f}')
    plt.text(ss[-1], aa[-1], f'{+qe[X2]:.4f}')
    
    #%% Dibujar los diagramas de fuerza cortante
    plt.figure(4)
    pos = T @ np.array([s      ,
                        esc_V*V]) # escalamiento del diagrama
    
    ss = pos[X,:] + x1
    vv = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-')
    if con_lineas:
        plt.plot(np.r_[x1, ss, x2], np.r_[y1, vv, y2], 'r-', linewidth=2)
    else:
        plt.plot(ss, vv, 'r-', linewidth=2)

    plt.text(ss[ 0], vv[ 0], f'{+qe[Y1]:.4f}')
    plt.text(ss[-1], vv[-1], f'{-qe[Y2]:.4f}')
    
    #%% Dibujar los diagramas de momento flector
    plt.figure(5)
    pos = T @ np.array([s      ,
                        esc_M*M]) # escalamiento del diagrama    
    
    ss = pos[X,:] + x1
    mm = pos[Y,:] + y1
    
    plt.plot([x1, x2], [y1, y2], 'b-')
    if con_lineas:
        plt.plot(np.r_[x1, ss, x2], np.r_[y1, mm, y2], 'r-', linewidth=2)
    else:
        plt.plot(ss, mm, 'r-', linewidth=2)
    plt.text(ss[ 0], mm[ 0], f'{-qe[M1]:.4f}')
    plt.text(ss[-1], mm[-1], f'{+qe[M2]:.4f}')

    # se grafican los máximos y los mínimos en caso que no estén en los bordes
    idminM = np.argmin(M); 
    if idminM not in [0, npuntos-1]: 
        plt.text(ss[idminM], mm[idminM], f'{M[idminM]:.4f}')
    
    idmaxM = np.argmax(M); 
    if idmaxM not in [0, npuntos-1]: 
        plt.text(ss[idmaxM], mm[idmaxM], f'{M[idmaxM]:.4f}')
