from leer_GMSH import xnod_from_msh, LaG_from_msh, plot_msh

malla = 'malla1.msh'

# Matriz de coordenadas nodales
xnod = xnod_from_msh(malla, dim=2)

# Matriz de interconexión nodal
mat_LaG = LaG_from_msh(malla)
mat = mat_LaG[:,0]
LaG = mat_LaG[:,1:]+1
nef = LaG.shape[0]

# Se grafica la malla:
plot_msh(malla, '2D', 
         mostrar_nodos=True, 
         mostrar_num_nodo=True, 
         mostrar_num_elem=True)


from obtener_grupos_fisicos import grupos_fisicos, obtener_nodos

# Obtener todos los grupos físicos de la malla:
dict_nombres, dict_nodos = grupos_fisicos(malla)

print('~'*50)
print('Grupos físicos reportados:\n')
for tag in dict_nombres.keys():
    dim    = dict_nombres[tag][0]
    nombre = dict_nombres[tag][1]
    nodos  = dict_nodos[tag]
    print(f'Grupo físico: {tag}\nDimensión: {dim}\n'
          f'Nombre: {nombre}\nContiene nodos: {nodos}.\n')

# Obtener nodos únicamente del grupo físico llamado "solido":
#nod_AB = obtener_nodos(malla, 'solido')
print('~'*50)
#print(f'Grupo físico "AB" contiene nodos: {nod_AB}')
