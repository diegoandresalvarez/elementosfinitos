## cargar variables desde archivos al sistema
# xnod - posicion de los nodos
xnod = readdlm("malla3_xnod.txt")
xnod = xnod[:,[2, 3]]
xnod = 0.01*xnod           # se pasa de cm a metros

# LaG - definicion de elementos finitos con respecto a nodos
LaG = readdlm("malla3_LaG.txt", Int)
LaG = LaG[:,2:9]                    # se borra la columna 1
LaG = LaG[:,vec([1 5 2 6 3 7 4 8])] # se renumeran los nodos


# Numeracion local GiD:    # Numeracion local Oñate:
#      ^ eta               #      ^ eta
#      |                   #      |
#      |                   #      |
#  4---7---3               #  7---6---5
#  |   |   |               #  |   |   |
#  8---+---6----> xi       #  8---+---4----> xi
#  |   |   |               #  |   |   |
#  1---5---2               #  1---2---3

## Se definen las restricciones
# X = 1; Y = 2;
nodos_restringidos = vec([1 3 5 9 12 18 24 33 43 50 66 80 99 123 151])
nnr = length(nodos_restringidos)
restricciones = [ nodos_restringidos fill(X,nnr) zeros(nnr)
                  nodos_restringidos fill(Y,nnr) zeros(nnr) ]
                  
