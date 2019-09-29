# presentación de resultados
import matplotlib.pyplot as plt

def imprimir_malla(nef, )
    plt.figure()
    cgx = np.zeros(nef); cgy = np.zeros(nef) # almacena el centro de gravedad
    for e in range(nef):
    idx = [NL1, NL2, NL3, NL1]
    plt.plot(xnod[LaG[e, idx], X], xnod[LaG[e, idx], Y], 'b')

    # Calculo la posición del centro de gravedad del triángulo
    cgx[e] = np.mean(xnod[LaG[e,:], X])
    cgy[e] = np.mean(xnod[LaG[e,:], Y])
    plt.text(cgx[e], cgy[e], f'{e+1}', horizontalalignment='center', 
                                        verticalalignment='center',   color='b')

    plt.plot(xnod[:,X], xnod[:,Y], 'r*')
    for i in range(nno):
    plt.text(xnod[i,X], xnod[i,Y], f'{i+1}', color='r')
    plt.axis('equal') # tight
    plt.title('Malla de elementos finitos')
    plt.show()