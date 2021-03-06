# Cálculo de una losa simplemente apoyada en sus cuatro bordes utilizando el EF de losa DKMQ (Discrete Kirchhoff Mindlin Quadrilateral)

La implementación del EF DKMQ está basada en los artículos:

> Katili, I. (1993), A new discrete Kirchhoff-Mindlin element based on Mindlin-Reissner plate theory and assumed shear strain fields-part II: An extended DKQ element for thick-plate bending analysis. Int. J. Numer. Meth. Engng., 36: 1885-1908. https://doi.org/10.1002/nme.1620361107

y

> Katili, et. al. (2018) - A comparative formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements with new numerical results based on s-norm tests. Computers & Structures, Volume 204, Pages 48-64, https://doi.org/10.1016/j.compstruc.2018.04.001.

Considere la losa mostrada en la figura
![](../../ejemplos/losa.png)

dicha losa tiene:
* dimensión: a = 2m, b = 4m, espesor t = 5 cm.
* material: E = 210 GPa y ν = 0.30.
* soporta una carga p = 10 kN/m^2, con u = 0.5m, v = 1m, ξ = 1.25m, η = 1.5 m.


La losa se calculó con elementos finitos DKMQ, implementado en el archivo [EF_DKMQ.m](EF_DKMQ.m), obteniendo los siguientes diagramas:

### Deformación vertical w

![figs/w.png](figs/w.png)

### Momentos flectores Mx y My y torsores Mxy

![figs/MxMyMxy.png](figs/MxMyMxy.png)

### Momentos flectores máximos y mínimos y momentos torsores máximos con sus respectivas inclinaciones

![figs/M1fM2fMtmax.png](figs/M1fM2fMtmax.png)

Recuerde que M1f y Mf2 nos dicen la forma de colocar el refuerzo óptimo en la parte superior e inferior de la losa, respectivamente.

### Momentos de diseño de Wood y Armer en la parte superior e inferior

![figs/WA_sup.png](figs/WA_sup.png)

![figs/WA_inf.png](figs/WA_inf.png)


### Comparación contra la solución teórica

El desplazamiento vertical de la losa se comparó contra su solución teórica, la cual se encuentra en Eduard Ventsel and Theodor Krauthammer (2001) - Thin plates and shells: theory, analysis and applications. Marcel Dekker : New York. Páginas 53 y 54.

La deformación teórica en la placa esta dada por:
<!---
Compile en: https://tex.s2cms.com

\begin{equation}
  w(x,y) = \frac{1}{\pi^4 D}\sum_{m=1}^\infty \sum_{n=1}^\infty
  \frac{p_{mn}}{\left(\frac{m^2}{a^2} + \frac{n^2}{b^2}\right)^2}
  \sin\left(\frac{m \pi x}{a}\right)
  \sin\left(\frac{n \pi y}{b}\right)
\end{equation}
donde
\begin{equation}
  p_{mn} = \frac{16 p}{\pi^2 m n}
  \sin\left(\frac{m \pi \xi}{a}\right)
  \sin\left(\frac{n \pi \eta}{b}\right)
  \sin\left(\frac{m \pi u}{2a}\right)
  \sin\left(\frac{n \pi v}{2b}\right)
\end{equation}
--->

![](https://i.upmath.me/svg/w(x%2Cy)%20%3D%20%5Cfrac%7B1%7D%7B%5Cpi%5E4%20D%7D%5Csum_%7Bm%3D1%7D%5E%5Cinfty%20%5Csum_%7Bn%3D1%7D%5E%5Cinfty%0A%20%20%5Cfrac%7Bp_%7Bmn%7D%7D%7B%5Cleft(%5Cfrac%7Bm%5E2%7D%7Ba%5E2%7D%20%2B%20%5Cfrac%7Bn%5E2%7D%7Bb%5E2%7D%5Cright)%5E2%7D%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bm%20%5Cpi%20x%7D%7Ba%7D%5Cright)%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bn%20%5Cpi%20y%7D%7Bb%7D%5Cright))

donde

![](https://i.upmath.me/svg/p_%7Bmn%7D%20%3D%20%5Cfrac%7B16%20p%7D%7B%5Cpi%5E2%20m%20n%7D%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bm%20%5Cpi%20%5Cxi%7D%7Ba%7D%5Cright)%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bn%20%5Cpi%20%5Ceta%7D%7Bb%7D%5Cright)%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bm%20%5Cpi%20u%7D%7B2a%7D%5Cright)%0A%20%20%5Csin%5Cleft(%5Cfrac%7Bn%20%5Cpi%20v%7D%7B2b%7D%5Cright))

obteniendo un error absoluto relativo máximo del 0.2128%. Las fórmulas anteriores se encuentran implementadas en el archivo [calc_w.m](calc_w.m).

Código elaborado por:
* Diego Andrés Alvarez Marín 
* Sebastián Jaramillo Moreno
