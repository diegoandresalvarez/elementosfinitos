# Ejemplo Uribe Escamilla 

# Unidades: toneladas, cm

# **** SE GENERA EL MODELO ****
# se borra la memoria
wipe                       

# se crea el modelo
# ndm -- numero de dimensiones
# ndf -- numero de grados de libertad por nodo
model BasicBuilder -ndm 2 -ndf 2

# se crean los nodos
# node nodeId coord_x coord_y
node 1   0.0   0.0
node 2 800.0   0.0
node 3 400.0 300.0
node 4 400.0   0.0

# se ponen las restricciones de los apoyos
# fix nodeID restring_x, restring_y (1=si, 0=no)
fix 1 1 1 
fix 2 0 1

# se definen las propiedades del material (se asume material elastico)
# uniaxialMaterial Elastic matID E
uniaxialMaterial Elastic 1 2040

# se definen los elementos de la cercha
# element truss trussID node1 node2 A matID
element Truss 1 1 3 100.0 1
element Truss 2 1 4  40.0 1 
element Truss 3 3 2 150.0 1
element Truss 4 4 2  40.0 1
element Truss 5 3 4  30.0 1
    

# **** SE DEFINEN LAS CARGAS ****    
# El factor de carga variara linealmente con el tiempo
# timeSeries Linear $tag
timeSeries Linear 1

set ang [expr atan2(3.0,4.0)]

# se utiliza un patron de cargas "Plain"
# pattern Plain $tag $timeSeriesTag { $loads }
pattern Plain 1 1 {    
    # se crean las cargas nodales
    # load nodeID xForce yForce
    load 3 [expr 5*cos($ang)] [expr 5*sin($ang)]
    load 4 0 -20    
}

# **** SE CONFIGURA COMO SE RESOLVERA EL SISTEMA Ka=f ****
# Create the system of equation, a SPD using a band storage scheme
system BandSPD

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the constraint handler, a Plain handler is used as homo constraints
constraints Plain

# Create the integration scheme, the LoadControl scheme using steps of 1.0
integrator LoadControl 1.0

# Create the solution algorithm, a Linear algorithm is created
algorithm Linear

# create the analysis object 
analysis Static 

# **** SE REALIZA EL ANALISIS ESTRUCTURAL ****
analyze 1

# **** SE IMPRIMEN LOS RESULTADOS ****
# se calculan las recciones
reactions

# se imprimen los desplazamientos del nodo 4
puts "*************************************************************************"
puts "Desplazamiento (x,y) del nodo 4: [nodeDisp 4]"

# se imprimen las reacciones del nodo 1
puts "*************************************************************************"
puts "Reacciones (Fx,Fy) del nodo 1: [nodeReaction 1]"

# se imprime la informacion sobre el nodo 4
puts "*************************************************************************"
print node 4

# se imprime la informacion sobre todas las barras (elementos)
puts "*************************************************************************"
print ele

# se imprime toda la informacion
puts "*************************************************************************"
print


# **** SE GUARDAN LOS RESULTADOS A DISCO ****
# create a Recorder object for the nodal displacements at node 4
recorder Node -file example.out -node 4 -dof 1 2 disp

# Create a recorder for element forces, one in global and the other local system
recorder Element -file eleGlobal.out -ele 1 2 3 forces
recorder Element -file eleLocal.out  -ele 1 2 3 basicForces

# se graban los resultados en un archivo
record