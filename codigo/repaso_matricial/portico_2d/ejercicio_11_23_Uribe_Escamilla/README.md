#  Ejemplo 11.23 de Uribe Escamilla: análisis de un pórtico bidimensional

Los programas:
* MATLAB: [ejemplo_portico.m](ejemplo_portico.m)
* PYTHON: [ejemplo_portico.py](ejemplo_portico.py)

estiman las reacciones en los apoyos, los desplazamientos y las fuerzas axiales de las barras del pórtico 2D que se muestra a continuación:

![figura](ejemplo_portico.svg)

Se incluye el archivo [ejemplo_portico.pdf](ejemplo_portico.pdf) para que se compare con el libro de Uribe Escamilla. Dicho libro se puede descargar libremente de:

https://www.researchgate.net/publication/31754481_Analisis_de_estructuras_J_Uribe_Escamilla (el ejemplo se encuentra en la página 526 de esta versión del libro).

De otro lado, los programas:
* MATLAB: [ejemplo_portico_resolv_ec_dif.m](ejemplo_portico_resolv_ec_dif.m)
* MATLAB: [fe_ec_dif.m](fe_ec_dif.m)
* MATLAB: [Ke_ec_dif.m](Ke_ec_dif.m)
resuelven el mismo problema pero crean la matriz de rigidez Ke y el vector de fuerzas nodales equivalentes fe mediante la solución de la ecuación diferencial.


