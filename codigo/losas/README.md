**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**
<span style="color: #0000ff;
font-size: 200%;">Nota: estos códigos están hechos para que funcionen con MATLAB 2013a. No los he actualizado a versiones más nuevas de MATLAB, ya que el 2013a es el MATLAB que se tiene instalado en los computadores de la Universidad Nacional de Colombia - Sede Manizales. Por lo tanto, los programas para álgebra simbólica podrían fallar si usted utiliza versiones modernas de MATLAB.</span>

[[image:http://imgs.xkcd.com/comics/ballmer_peak.png]]
Fuente: [[http://xkcd.com/323/]]



=CAPITULO 8=
==Programa para obtener las funciones de forma del elemento finito de losa MZC==
** Código MATLAB: [[file:c8_func_forma_MZC.m]] **usa el toolbox de álgebra simbolica**
** Código compatible con MATLAB 2013a: [[file:c8_func_forma_MZC_MATLAB_2013a.m]] **usa el toolbox de álgebra simbolica**

Con este código se obtuvieron, por ejemplo, las siguientes funciones de forma:
[[code]]
                                  2           2
           (eta - 1) (xi + 1) (eta  + eta + xi  - xi - 2)
N{2}   =   ----------------------------------------------
                                  8


             (eta - 1) (xi - 1) (xi + 1)
Nb{2}  =   - ----------------------------
                          8        

                   2
          (eta - 1)  (eta + 1) (xi + 1)
Nbb{2} =  -----------------------------
                        8
[[code]]

Siendo los gráficos de estas funciones:
[[image:c8_func_forma_MZC.png]]

Con el mismo programa se obtiene la matriz de rigidez del elemento (una matriz de 12x12) que se muestra a continuación
[[math]]
\Large
\boldsymbol{K}^{(e)} = \frac{D}{ab}
\left(\begin{smallmatrix}
\frac{b^2}{a^2} - \frac{\nu}{5} + \frac{a^2}{b^2} + \frac{7}{10} & \frac{2 \nu}{5} + \frac{b^2}{a^2} + \frac{1}{10} & \frac{2 \nu}{5} + \frac{a^2}{b^2} + \frac{1}{10} & \frac{\nu}{5} - \frac{b^2}{a^2} + \frac{a^2}{2 b^2} - \frac{7}{10} & \frac{b^2}{a^2} - \frac{\nu}{10} + \frac{1}{10} & \frac{a^2}{2 b^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{7}{10} - \frac{b^2}{2 a^2} - \frac{a^2}{2 b^2} - \frac{\nu}{5} & \frac{\nu}{10} + \frac{b^2}{2 a^2} - \frac{1}{10} & \frac{\nu}{10} + \frac{a^2}{2 b^2} - \frac{1}{10} & \frac{\nu}{5} + \frac{b^2}{2 a^2} - \frac{a^2}{b^2} - \frac{7}{10} & \frac{b^2}{2 a^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{a^2}{b^2} - \frac{\nu}{10} + \frac{1}{10}\\
\frac{2 \nu}{5} + \frac{b^2}{a^2} + \frac{1}{10} & \frac{4 b^2}{3 a^2} - \frac{4 \nu}{15} + \frac{4}{15} & \nu & \frac{\nu}{10} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{\nu}{15} + \frac{2 b^2}{3 a^2} - \frac{1}{15} & 0 & \frac{1}{10} - \frac{b^2}{2 a^2} - \frac{\nu}{10} & \frac{b^2}{3 a^2} - \frac{\nu}{15} + \frac{1}{15} & 0 & \frac{b^2}{2 a^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{4 \nu}{15} + \frac{2 b^2}{3 a^2} - \frac{4}{15} & 0\\
\frac{2 \nu}{5} + \frac{a^2}{b^2} + \frac{1}{10} & \nu & \frac{4 a^2}{3 b^2} - \frac{4 \nu}{15} + \frac{4}{15} & \frac{a^2}{2 b^2} - \frac{2 \nu}{5} - \frac{1}{10} & 0 & \frac{4 \nu}{15} + \frac{2 a^2}{3 b^2} - \frac{4}{15} & \frac{1}{10} - \frac{a^2}{2 b^2} - \frac{\nu}{10} & 0 & \frac{a^2}{3 b^2} - \frac{\nu}{15} + \frac{1}{15} & \frac{\nu}{10} - \frac{a^2}{b^2} - \frac{1}{10} & 0 & \frac{\nu}{15} + \frac{2 a^2}{3 b^2} - \frac{1}{15}\\
\frac{\nu}{5} - \frac{b^2}{a^2} + \frac{a^2}{2 b^2} - \frac{7}{10} & \frac{\nu}{10} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{a^2}{2 b^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{b^2}{a^2} - \frac{\nu}{5} + \frac{a^2}{b^2} + \frac{7}{10} &  - \frac{2 \nu}{5} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{2 \nu}{5} + \frac{a^2}{b^2} + \frac{1}{10} & \frac{\nu}{5} + \frac{b^2}{2 a^2} - \frac{a^2}{b^2} - \frac{7}{10} & \frac{2 \nu}{5} - \frac{b^2}{2 a^2} + \frac{1}{10} & \frac{a^2}{b^2} - \frac{\nu}{10} + \frac{1}{10} & \frac{7}{10} - \frac{b^2}{2 a^2} - \frac{a^2}{2 b^2} - \frac{\nu}{5} & \frac{1}{10} - \frac{b^2}{2 a^2} - \frac{\nu}{10} & \frac{\nu}{10} + \frac{a^2}{2 b^2} - \frac{1}{10}\\
\frac{b^2}{a^2} - \frac{\nu}{10} + \frac{1}{10} & \frac{\nu}{15} + \frac{2 b^2}{3 a^2} - \frac{1}{15} & 0 &  - \frac{2 \nu}{5} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{4 b^2}{3 a^2} - \frac{4 \nu}{15} + \frac{4}{15} & - \nu & \frac{2 \nu}{5} - \frac{b^2}{2 a^2} + \frac{1}{10} & \frac{4 \nu}{15} + \frac{2 b^2}{3 a^2} - \frac{4}{15} & 0 & \frac{\nu}{10} + \frac{b^2}{2 a^2} - \frac{1}{10} & \frac{b^2}{3 a^2} - \frac{\nu}{15} + \frac{1}{15} & 0\\
\frac{a^2}{2 b^2} - \frac{2 \nu}{5} - \frac{1}{10} & 0 & \frac{4 \nu}{15} + \frac{2 a^2}{3 b^2} - \frac{4}{15} & \frac{2 \nu}{5} + \frac{a^2}{b^2} + \frac{1}{10} & - \nu & \frac{4 a^2}{3 b^2} - \frac{4 \nu}{15} + \frac{4}{15} & \frac{\nu}{10} - \frac{a^2}{b^2} - \frac{1}{10} & 0 & \frac{\nu}{15} + \frac{2 a^2}{3 b^2} - \frac{1}{15} & \frac{1}{10} - \frac{a^2}{2 b^2} - \frac{\nu}{10} & 0 & \frac{a^2}{3 b^2} - \frac{\nu}{15} + \frac{1}{15}\\

\frac{7}{10} - \frac{b^2}{2 a^2} - \frac{a^2}{2 b^2} - \frac{\nu}{5} & \frac{1}{10} - \frac{b^2}{2 a^2} - \frac{\nu}{10} & \frac{1}{10} - \frac{a^2}{2 b^2} - \frac{\nu}{10} & \frac{\nu}{5} + \frac{b^2}{2 a^2} - \frac{a^2}{b^2} - \frac{7}{10} & \frac{2 \nu}{5} - \frac{b^2}{2 a^2} + \frac{1}{10} & \frac{\nu}{10} - \frac{a^2}{b^2} - \frac{1}{10} & \frac{b^2}{a^2} - \frac{\nu}{5} + \frac{a^2}{b^2} + \frac{7}{10} &  - \frac{2 \nu}{5} - \frac{b^2}{a^2} - \frac{1}{10} &  - \frac{2 \nu}{5} - \frac{a^2}{b^2} - \frac{1}{10} & \frac{\nu}{5} - \frac{b^2}{a^2} + \frac{a^2}{2 b^2} - \frac{7}{10} & \frac{\nu}{10} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{2 \nu}{5} - \frac{a^2}{2 b^2} + \frac{1}{10}\\
\frac{\nu}{10} + \frac{b^2}{2 a^2} - \frac{1}{10} & \frac{b^2}{3 a^2} - \frac{\nu}{15} + \frac{1}{15} & 0 & \frac{2 \nu}{5} - \frac{b^2}{2 a^2} + \frac{1}{10} & \frac{4 \nu}{15} + \frac{2 b^2}{3 a^2} - \frac{4}{15} & 0 &  - \frac{2 \nu}{5} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{4 b^2}{3 a^2} - \frac{4 \nu}{15} + \frac{4}{15} & \nu & \frac{b^2}{a^2} - \frac{\nu}{10} + \frac{1}{10} & \frac{\nu}{15} + \frac{2 b^2}{3 a^2} - \frac{1}{15} & 0\\
 \frac{\nu}{10} + \frac{a^2}{2 b^2} - \frac{1}{10} & 0 & \frac{a^2}{3 b^2} - \frac{\nu}{15} + \frac{1}{15} & \frac{a^2}{b^2} - \frac{\nu}{10} + \frac{1}{10} & 0 & \frac{\nu}{15} + \frac{2 a^2}{3 b^2} - \frac{1}{15} &  - \frac{2 \nu}{5} - \frac{a^2}{b^2} - \frac{1}{10} & \nu & \frac{4 a^2}{3 b^2} - \frac{4 \nu}{15} + \frac{4}{15} & \frac{2 \nu}{5} - \frac{a^2}{2 b^2} + \frac{1}{10} & 0 & \frac{4 \nu}{15} + \frac{2 a^2}{3 b^2} - \frac{4}{15}\\
 \frac{\nu}{5} + \frac{b^2}{2 a^2} - \frac{a^2}{b^2} - \frac{7}{10} & \frac{b^2}{2 a^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{\nu}{10} - \frac{a^2}{b^2} - \frac{1}{10} & \frac{7}{10} - \frac{b^2}{2 a^2} - \frac{a^2}{2 b^2} - \frac{\nu}{5} & \frac{\nu}{10} + \frac{b^2}{2 a^2} - \frac{1}{10} & \frac{1}{10} - \frac{a^2}{2 b^2} - \frac{\nu}{10} & \frac{\nu}{5} - \frac{b^2}{a^2} + \frac{a^2}{2 b^2} - \frac{7}{10} & \frac{b^2}{a^2} - \frac{\nu}{10} + \frac{1}{10} & \frac{2 \nu}{5} - \frac{a^2}{2 b^2} + \frac{1}{10} & \frac{b^2}{a^2} - \frac{\nu}{5} + \frac{a^2}{b^2} + \frac{7}{10} & \frac{2 \nu}{5} + \frac{b^2}{a^2} + \frac{1}{10} &  - \frac{2 \nu}{5} - \frac{a^2}{b^2} - \frac{1}{10}\\
 \frac{b^2}{2 a^2} - \frac{2 \nu}{5} - \frac{1}{10} & \frac{4 \nu}{15} + \frac{2 b^2}{3 a^2} - \frac{4}{15} & 0 & \frac{1}{10} - \frac{b^2}{2 a^2} - \frac{\nu}{10} & \frac{b^2}{3 a^2} - \frac{\nu}{15} + \frac{1}{15} & 0 & \frac{\nu}{10} - \frac{b^2}{a^2} - \frac{1}{10} & \frac{\nu}{15} + \frac{2 b^2}{3 a^2} - \frac{1}{15} & 0 & \frac{2 \nu}{5} + \frac{b^2}{a^2} + \frac{1}{10} & \frac{4 b^2}{3 a^2} - \frac{4 \nu}{15} + \frac{4}{15} & - \nu\\
 \frac{a^2}{b^2} - \frac{\nu}{10} + \frac{1}{10} & 0 & \frac{\nu}{15} + \frac{2 a^2}{3 b^2} - \frac{1}{15} & \frac{\nu}{10} + \frac{a^2}{2 b^2} - \frac{1}{10} & 0 & \frac{a^2}{3 b^2} - \frac{\nu}{15} + \frac{1}{15} & \frac{2 \nu}{5} - \frac{a^2}{2 b^2} + \frac{1}{10} & 0 & \frac{4 \nu}{15} + \frac{2 a^2}{3 b^2} - \frac{4}{15} &  - \frac{2 \nu}{5} - \frac{a^2}{b^2} - \frac{1}{10} & - \nu & \frac{4 a^2}{3 b^2} - \frac{4 \nu}{15} + \frac{4}{15}
\end{smallmatrix}\right)
[[math]]

donde

[[math]]
D = \frac{E t^3}{12(1-\nu)}
[[math]]

==Ejemplo de cálculo de la deflexión en una losa rectangular utilizando elementos finitos MZC==
Considere la losa mostrada en la Figura
[[image:c8_losa_rect.png]]
donde 
[[math]]
\\
a = 2 \text{ m}\\
b = 4 \text{ m}\\
\xi = 1.25 \text{ m}\\
\eta = 1.5 \text{ m}\\
u = 0.5 \text{ m}\\
v = 1 \text{ m}\\
t = 0.05 \text{ m}\\
t  = 0.05 \text{ espesor de la losa (m)}\\
p  = -10 \text{ kN/m}^2 \text{ (carga)}\\
E = 210 \text{ GPa (m\'odulo de Young)} \\
\nu = 0.3 \text{ (coeficiente de Poisson)}
[[math]]

Con el programa [[file:c8_losa_rect_MZC.zip]] se calcularon las deflexiones y los momentos en la losa, como se muestra a continuación:

Las deflexiones son:
[[image:c8_deformada_losa_rect.png]]

Los momentos son: <span style="color: #ff0000;">(ESTOS SE CALCULARON EN EL CENTRO DEL ELEMENTO, SIN EMBARGO ES MEJOR CALCULARLOS EN LOS PUNTOS DE GAUSS 2x2 (VER OÑATE, PAG 341). FALTA HACERLO)</span>
[[image:c9_momentos_losa_rect.png]]

La solución analítica a este problema está dado por:
[[math]]
   w(x,y) = \frac{1}{\pi^4 D}\sum_{m=1}^\infty \sum_{n=1}^\infty
   \frac{p_{mn}}{\left(\frac{m^2}{a^2} + \frac{n^2}{b^2}\right)^2}
   \sin\left(\frac{m \pi x}{a}\right)
   \sin\left(\frac{n \pi y}{b}\right)
[[math]]
donde
[[math]]
   p_{mn} = \frac{16 p}{\pi^2 m n}
   \sin\left(\frac{m \pi \xi}{a}\right)
   \sin\left(\frac{n \pi \eta}{b}\right)
   \sin\left(\frac{m \pi u}{2a}\right)
   \sin\left(\frac{n \pi v}{2b}\right)
[[math]]

En el codigo anterior se programó dicha ecuación y se observa la precisión del método de los elementos finitos comparado con la respuesta analíticamente exacta. De hecho se observó que para la malla utilizada el error relativo está máximo es del 0.28%; es decir, los errores son extremadamente pequeños!!!


==Calculo de las funciones de forma y matriz L del elemento de Tocher==
El programa [[file:c8_func_forma_Tocher.m]] deduce y grafica las mencionadas funciones de forma. Por ejemplo, las funciones de forma del nodo 3 son:
[[image:c8_func_forma_Tocher.png]]


**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**
<span style="color: #0000ff;
font-size: 200%;">Nota: estos códigos están hechos para que funcionen con MATLAB 2013a. No los he actualizado a versiones más nuevas de MATLAB, ya que el 2013a es el MATLAB que se tiene instalado en los computadores de la Universidad Nacional de Colombia - Sede Manizales. Por lo tanto, los programas para álgebra simbólica podrían fallar si usted utiliza versiones modernas de MATLAB.</span>

[[image:http://imgs.xkcd.com/comics/ballmer_peak.png]]
Fuente: [[http://xkcd.com/323/]]

=CAPITULO 9: Placas gruesas: teoría de Reissner-Mindlin=

==Elemento isoparametrico QL9 (solucion a una losa y ejemplo con modo de energía nula propagable -- Integración reducida)==
Código de MATLAB: [[file:c9_QL9_integracion_reducida.zip]]
[[image:c9_MEN_QL9.png]]

==Elemento QHET (elemento heterosis)==
Código de MATLAB: [[file:c9_QHET_elemento_heterosis.zip]]

==Ejemplo Oñate p. 320 (inglés), p. 388 (español): imposición de las deformaciones angulares==
Código de MATLAB: [[file:c9_motivacion_ANS.m]]

==Calculo de las matrices A*inv(P)*T para los elementos finitos QLLL, QQQQ-L, QQQQ-S y TQQL (imposición del campo vectorial de deformaciones angulares)==
Código de MATLAB: [[file:c9_APm1T_QLLL___QQQQ_S___QQQQ_L___TQQL.zip]]

==Elemento de losa de RM QQQQ_L (imposición del campo vectorial de deformaciones angulares)==
Código de MATLAB: [[file:c9_QQQQ_L_Imposicion_campo_deformaciones_cortantes.zip]]
