##  Análisis de 3 barras trabajando a tracción

Los programas de esta carpeta, estiman las reacciones en los apoyos (nodos 1 y 2) y los desplazamientos en los nodos no restringidos (nodos 3 y 4) de la figura mostrada:

![figura](tres_barras_a_traccion.svg)

Así mismo, `ejemplo_barra_v2` estima las fuerzas axiales en cada una de las barras.

Los programas de MATLAB (*.m) resuelven el problema utilizando el **toolbox de álgebra simbólica**.

La versión de PYTHON (*.py) utiliza **el módulo SymPy**.

* `ejemplo_barra_v1`: se le especifica "a mano" la matriz *K* y resuelve el sistema de ecuaciones *Ka - f = q*.
* `ejemplo_barra_v2`: este programa ensambla la matriz *K* y resuelve el sistema de ecuaciones *Ka - f = q*.

La solución analítica de este ejemplo se encuentra [aquí](../../../diapositivas/01_demo_q_Ka-f_en_barra.pdf).
