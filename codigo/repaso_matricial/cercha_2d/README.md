#  Ejemplo 11.3 Uribe Escamilla: análisis de una cercha

Los programas de esta carpeta, estiman las reacciones en los apoyos, los desplazamientos y las fuerzas axiales de las barras de las cerchas que se muestran a continuación:

## Cercha sin apoyo inclinado
![figura](ej_11_3_uribe_escamilla.svg)

El material es acero, el cual tiene un módulo de elasticidad *E* de 2040 ton/cm². Las áreas *A* están entre paréntesis en cm².

Se incluye el archivo `.pdf` para que se compare con el libro de Uribe Escamilla. Dicho libro se puede descargar libremente de:

https://www.researchgate.net/publication/31754481_Analisis_de_estructuras_J_Uribe_Escamilla (el ejemplo se encuentra en la página 437 de esta versión del libro).

Souciones con:
* MATLAB: [ejemplo_cercha.m](ejemplo_cercha.m)
* PYTHON: [ejemplo_cercha.py](ejemplo_cercha.py)
* [OpenSees](https://opensees.berkeley.edu/): [ejemplo_cercha_OpenSees.tcl](ejemplo_cercha_OpenSees.tcl)

## Cercha con apoyo inclinado
![figura](ejemplo_cercha_apoyo_inclinado.svg)

El material es acero, el cual tiene un módulo de elasticidad *E* de 2040 ton/cm². Las áreas *A* están entre paréntesis en cm².

Solución con:
* MATLAB: [ejemplo_cercha_apoyo_inclinado.m](ejemplo_cercha_apoyo_inclinado.m)
