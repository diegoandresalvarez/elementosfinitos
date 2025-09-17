import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# %%
def timoshenko_ode(x, y, L, EI, GA_star, w1, w2):
    """
    Define el sistema de Ecuaciones Diferenciales Ordinarias (EDO) para la viga
    de Timoshenko.
    
    El vector de estado y es:
    y[0] = V   (fuerza cortante)
    y[1] = M   (momento flector)
    y[2] = θ   (rotación)
    y[3] = w   (deflexión)
    
    El sistema de EDOs es y' = f(x, y):
    dV/dx = q(x)
    dM/dx = V
    dθ/dx = M/EI
    dw/dx = θ - V/GA_star
    """
    # Indices para mayor la legibilidad
    V, M, θ, w = 0, 1, 2, 3

    # Carga trapezoidal q(x)
    q = w1 + (w2 - w1)*x/L
    
    # Derivadas
    dV_dx = q
    dM_dx = y[V]
    dθ_dx = y[M]/EI
    dw_dx = y[θ] - y[V]/GA_star

    return np.vstack([dV_dx, dM_dx, dθ_dx, dw_dx])

# %%
def calcular_matriz_rigidez_num(L, EI, GA_star):
    """
    Calcula la matriz de rigidez 4x4 para el elemento de viga de Timoshenko.
    
    Lo hace resolviendo el BVP para 4 casos de condiciones de borde,
    correspondientes a desplazamientos/giros unitarios en cada grado de libertad.
    """
    # Indices para mayor la legibilidad
    V, M, θ, w = 0, 1, 2, 3

    K = np.zeros((4, 4))
    
    # Malla de puntos para el solucionador (de 0 a L)
    x_mesh = np.linspace(0, L, 10)
    
    # Bucle sobre las 4 condiciones de borde (4 columnas de la matriz K)
    for i in range(4):
        # Función para las condiciones de borde (Boundary Conditions - BC)
        def bc_gdl_i_en_1(y0, yL):
            # y0: solución en x=0,  yL: solución en x=L
            # Condiciones de borde base (empotrado-empotrado)
            res = np.array([ y0[w] - 0,      # w(0)=0
                             y0[θ] - 0,      # θ(0)=0
                             yL[w] - 0,      # w(L)=0
                             yL[θ] - 0 ])    # θ(L)=0            

            # Imponer desplazamiento/giro unitario en el grado de libertad 'i'
            res[i] = res[i] - 1.0
            return res

        # Estimación inicial de la solución (un vector de ceros es suficiente)
        y_guess = np.zeros((4, x_mesh.size))

        # La función de la EDO se pasa con lambda para fijar los parámetros y se
        # impone el caso sin carga w1=0, w2=0
        ode = lambda x, y: timoshenko_ode(x, y, L, EI, GA_star, w1=0, w2=0)

        # Resolver el BVP para el caso sin carga (w1=0, w2=0)
        sol = solve_bvp(ode, bc_gdl_i_en_1, x_mesh, y_guess)
        
        # Extraer las reacciones de la solución en los bordes
        V0 = sol.y[V,  0]  # Cortante en x=0
        M0 = sol.y[M,  0]  # Momento  en x=0
        VL = sol.y[V, -1]  # Cortante en x=L
        ML = sol.y[M, -1]  # Momento  en x=L

        # Ensamblar la columna 'i' de la matriz de rigidez
        K[:, i] = [+V0, -M0, -VL, +ML]
        
    return K

# %%
def calcular_vector_fuerzas_nodales_equiv_num(L, EI, GA_star, w1, w2):
    """
    Calcula el vector de fuerzas nodales equivalentes para una carga trapezoidal.
    
    Resuelve el BVP para una viga empotrada-empotrada bajo la carga dada.
    """
    # Indices para mayor la legibilidad
    V, M, θ, w = 0, 1, 2, 3

    # Malla de puntos
    x_mesh = np.linspace(0, L, 20)
    
    # Condiciones de borde para una viga empotrada-empotrada:
    def bc_empotrada_empotrada(y0, yL):
        return np.array([ y0[w] - 0,      # w(0)=0
                          y0[θ] - 0,      # θ(0)=0
                          yL[w] - 0,      # w(L)=0
                          yL[θ] - 0 ])    # θ(L)=0

    # Estimación inicial de la solución
    y_guess = np.zeros((4, x_mesh.size))

    # La función de la EDO se pasa con lambda para fijar los parámetros
    ode = lambda x, y: timoshenko_ode(x, y, L, EI, GA_star, w1, w2)
    
    # Resolver el BVP con la carga trapezoidal
    sol = solve_bvp(ode, bc_empotrada_empotrada, x_mesh, y_guess)

    # Extraer las reacciones en los apoyos
    V0 = sol.y[V,  0]  # Cortante en x=0
    M0 = sol.y[M,  0]  # Momento  en x=0
    VL = sol.y[V, -1]  # Cortante en x=L
    ML = sol.y[M, -1]  # Momento  en x=L    
    
    # El vector de fuerzas nodales equivalentes es el negativo de las reacciones
    # f1y = -V(0), f1m = +M(0), f2y = -(-V(L)), f2m = -(+M(L))
    f = np.array([-V0, +M0, +VL, -ML])
    
    return f, sol

# %%
def calcular_matriz_rigidez_exacta(L, EI, GA_star):
    """
    Calcula la matriz de rigidez exacta para el elemento de viga de Timoshenko.
    """
    # Parámetro beta
    beta = (12*EI)/(L**2 * GA_star)
    
    # Matriz de rigidez de flexión del elemento
    Ke = (EI/((1 + beta)*L**3)) * np.array([
        [12,                6*L,    -12,               6*L],
        [6*L,   (4 + beta)*L**2,   -6*L,   (2 - beta)*L**2],
        [-12,              -6*L,     12,              -6*L],
        [6*L,   (2 - beta)*L**2,   -6*L,   (4 + beta)*L**2]
    ])
    
    return Ke

# %%    
def calcular_vector_fuerzas_exacto(L, EI, GA_star, w1, w2):
    """
    Calcula el vector de fuerzas nodales equivalentes exacto para una carga trapezoidal.
    """
    # Vector de fuerzas nodales equivalentes de una carga trapezoidal
    fe = np.array([
        # Y1
        (L * (80*EI*w1 + 40*EI*w2 + 7*GA_star*L**2*w1 + 3*GA_star*L**2*w2)) / 
        (20*GA_star*L**2 + 240*EI),
        
        # M1  
        (L**2 * (30*EI*w1 + 30*EI*w2 + 3*GA_star*L**2*w1 + 2*GA_star*L**2*w2)) /
        (60*(GA_star*L**2 + 12*EI)),
        
        # Y2
        (L * (40*EI*w1 + 80*EI*w2 + 3*GA_star*L**2*w1 + 7*GA_star*L**2*w2)) /
        (20*(GA_star*L**2 + 12*EI)),
        
        # M2
        -(L**2 * (30*EI*w1 + 30*EI*w2 + 2*GA_star*L**2*w1 + 3*GA_star*L**2*w2)) /
        (60*(GA_star*L**2 + 12*EI))
    ])
    
    return fe

# %% --- PROGRAMA PRINCIPAL ---
if __name__ == '__main__':
    # Carga trapezoidal
    w1 = -10000              # [N/m] Carga distribuida en x=0
    w2 = -30000              # [N/m] Carga distribuida en x=L

    # Viga de acero
    E  = 200e9               # [Pa]  Módulo de Young
    nu = 0.3                 # [-]   Coeficiente de Poisson
    G = E/(2*(1 + nu))       # [Pa]  Módulo de Cortante
    
    # Viga rectangular
    L = 5.0                  # [m]   Longitud de la viga
    h = 0.300                # [m]   Altura
    b = 0.150                # [m]   Ancho

    # Propiedades de la sección
    I       = b*h**3/12      # [m^4]    Momento de inercia
    A       = b*h            # [m^2]    Área
    kappa   = 5/6            # [-]      Factor de corrección de cortante
    A_star  = kappa*A        # [m^2]    Área efectiva a cortante A*
    EI      = E*I            # [Pa*m^4] Rigidez a flexión
    GA_star = G*A_star       # [Pa*m^2] Rigidez a cortante

    # Calcular la matriz de rigidez y el vector de fuerzas nodales equivalentes
    # de forma numérica
    K_num      = calcular_matriz_rigidez_num(L, EI, GA_star)
    f_num, sol = calcular_vector_fuerzas_nodales_equiv_num(L, EI, GA_star, w1, w2)

    # con las fórmulas exactas 
    K_exacta = calcular_matriz_rigidez_exacta(L, EI, GA_star)
    f_exacto = calcular_vector_fuerzas_exacto(L, EI, GA_star, w1, w2)

    # Imprimir los resultados
    np.set_printoptions(linewidth=200)

    # Matriz de rigidez
    print("\nMatriz de Rigidez K_TE (Numérica):")
    print(K_num)
    
    print("\nMatriz de Rigidez K_TE (Exacta):")
    print(K_exacta)

    print("\nDiferencia K_exacta - K_numérica:")
    print(K_exacta - K_num)

    K_iguales = np.allclose(K_exacta, K_num)
    print("\nSe pueden considerar las matrices K iguales:", K_iguales)
    if not K_iguales:
        print("Incremente el número de puntos en x_mesh para mejorar la precisión.")
    
    # Vector de fuerzas nodales equivalentes
    print("\nVector de Fuerzas Nodales Equivalentes f_TE (Numérico):")
    print(f_num)
    
    print("\nVector de Fuerzas Nodales Equivalentes f_TE (Exacto):")
    print(f_exacto)
    
    print("\nDiferencia f_exacto - f_numérico:")
    print(f_exacto - f_num)

    f_iguales = np.allclose(f_exacto, f_num)
    print("\nSe pueden considerar los vectores f iguales:", f_iguales)
    if not f_iguales:
        print("Incremente el número de puntos en x_mesh para mejorar la precisión.")
    
    # Graficar la deflexión y el momento flector bajo la carga
    x_plot = np.linspace(0, L, 200)
    y_plot = sol.sol(x_plot) # Evaluar la solución en más puntos
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Gráfico de Deflexión
    ax1.plot(x_plot, 1000*y_plot[0], label='Deflexión w(x)', color='b')
    ax1.set_ylabel('Deflexión (mm)')
    ax1.set_title(f'Análisis de Viga de Timoshenko (L={L}m, carga {w1/1e3}-{w2/1e3} kN/m)')
    ax1.grid(True)
    
    # Gráfico de Momento Flector
    ax2.plot(x_plot, y_plot[2]/1000, label='Momento M(x)', color='r')
    ax2.set_xlabel('Posición x (m)')
    ax2.set_ylabel('Momento Flector (kN·m)')
    ax2.grid(True)
    
    plt.tight_layout()
    plt.show()