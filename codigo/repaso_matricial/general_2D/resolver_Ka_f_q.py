def crear_C_g_enlace_rigido(nodo_padre, nodo_hijo):
    
    
    return C, g


def resolver_Ka_f_q(K, f, c, ac, C=None, g=None, idx_gdl_hijos=None):
    '''
    Esta función resuelve el sistema de ecuaciones Ka - f = q

    PARAMETROS DE ENTRADA:
    K  = matriz de rigidez global
    f  = vector de fuerzas nodales equivalentes    
    c  = GDL del desplazamiento conocidos (GDL de los apoyos)
    ac = desplazamientos conocidos

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
        # grados de libertad del desplazamiento desconocidos
        d = np.setdiff1d(np.arange(ngdl), c)
        
        # %% si no hay restricciones/enlaces rígidos:
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
        # %%si hay restricciones/enlaces rígidos
        # METODO 1: transformacion de K*a-f = q
        e = idx_gdl_hijos # idx de los GDL asociados a los nodos hijos
        
        # se calculan los idx de los GDL asociados a los nodos a retener
        # (incluye los nodos maestros)
        r = np.setdiff1d(np.arange(ngdl), e)
        
        # se verifica que ningún GDL de desplazamiento conocido sea un GDL hijo
        if np.any(np.intersect1d(c, e)):
            raise('No pueden haber GDL del empotramiento en los GDL hijos')
        
        num_hijos = C.shape[0]         # número de GDL hijos
        n_c       = ngdl - num_hijos   # número de GDL a retener
       
        Cr = C[:,r];   Ce = C[:,e]
        
        # Se verifica si Ce es singular o no
        if np.abs(np.linalg.det(Ce)) < 1e-5:
            raise('Ce debe ser invertible')
        
        I  = np.eye(n_c)
        T  = np.r_[  I,
                    -np.linalg.solve(Ce, Cr) ] # = -inv(Ce)*Cr
        
        O  = np.zeros(n_c)
        g0 = np.r_[  O,
                     np.linalg.solve(Ce, g)  ] # =  inv(Ce)*g
        
        # se reordenan los GDL de K y de f
        idx_re = np.c_[r, e]
        
        Kre = K[np.ix_(idx_re, idx_re)]
        fre = f[idx_re]
           
        Kr = T.T @ Kre @ T
        fr = T.T @ (fre - Kre@g0)
        
        # Se definen GDL conocidos y desconocidos asociados a los desplazamientos
        # idx con la numeración original (del grafico)
        # observe que en Kr*ar - fr = qr no hay nodos esclavos
        d = np.setdiff1d(r, c)
        
        # idx con la numeracion del reordenamiento [r e], esto es idx_re
        c_re = [ r.index(i) for i in c ]
        d_re = [ r.index(i) for i in d ]
        
        # Se descomponen los vectores a, f y la matriz K
        Krcc = Kr[np.ix_(c_re,c_re)];      Krcd = Kr[np.ix_(c_re,d_re)]
        Krdc = Kr[np.ix_(d_re,c_re)];      Krdd = Kr[np.ix_(d_re,d_re)]
        frc  = fr[d_re];                   frd  = fr[c_re]
        
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
        a_re = T@ar + g0
        
        # se calcula el vector q_re
        q_re = np.zeros(ngdl)  # fuerzas nodales equivalentes
        q_re[c_re] = qrd
        
        # se reordena a de la indexación re a la indexación original (del grafico)
        idx_orig  = [ idx_re.index(i) for i in range(ngdl) ]
        a = a_re[idx_orig]
        q = q_re[idx_orig]
        
    # %%
    return a, q
