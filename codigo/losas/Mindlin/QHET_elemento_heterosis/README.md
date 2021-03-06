# Teoría de losas Mindlin. Elemento finito heterosis QHET

## Cálculo de una losa simplemente apoyada en sus cuatro bordes utilizando el EF de losa heterosis QHET

Considere la losa mostrada en la figura

![](../../ejemplos/losa.png)

dicha losa tiene:
* dimensión: a = 2m, b = 4m, espesor t = 5 cm.
* material: E = 210 GPa y ν = 0.30.
* soporta una carga p = 10 kN/m^2, con u = 0.5m, v = 1m, ξ = 1.25m, η = 1.5 m.


La losa se calculó con elementos finitos QL9 en el archivo [ejemplo_losa_heterosis_QHET.m](ejemplo_losa_heterosis_QHET.m) obteniendo los siguientes diagramas:

### Deformación vertical w

![figs/w.png](figs/w.png)

### Momentos flectores Mx y My y torsores Mxy

![figs/MxMyMxy.png](figs/MxMyMxy.png)

### Momentos flectores máximos y mínimos y momentos torsores máximos con sus respectivas inclinaciones

![figs/M1fM2fMtmax.png](figs/M1fM2fMtmax.png)

Recuerde que M1f y Mf2 nos dicen la forma de colocar el refuerzo óptimo en la parte superior e inferior de la losa, respectivamente.

### Fuerzas cortantes Qx y Qy

![figs/QxQy.png](figs/QxQy.png)

### Fuerzas cortantes máximas Qmax con su respectivas inclinaciones

![figs/Qmax.png](figs/Qmax.png)
