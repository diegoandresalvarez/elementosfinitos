#  Análisis de un pórticos con rótulas

El programa [deduccion_Ke_empotrado_rodillo.m](deduccion_Ke_empotrado_rodillo.m) sirve para calcular las matrices de rigidez de elementos de pórtico que tienen a uno de sus lados una rótula.

Por ejemplo, la matriz de rigidez para:

![](figs/Ke_empotrado_rotula.svg)

es:
<!---
Compile en: https://tex.s2cms.com

\renewcommand\arraystretch{1.4}
\begin{bmatrix}
X_i\\
Y_i\\
M_i\\
X_j\\
Y_j\\
\end{bmatrix}
=
\begin{bmatrix}
\frac{EA}{L} &   0               &  0                & -\frac{EA}{L} &  0               \\
 0            &  \frac{3 EI}{L^3} & \frac{3 EI}{L^2}  &  0            & -\frac{3 EI}{L^3}\\
 0            &  \frac{3 EI}{L^2} & \frac{3 EI}{L}    &  0            & -\frac{3 EI}{L^2}\\
-\frac{EA}{L} &  0                &  0                & \frac{EA}{L}  &  0               \\
 0            & -\frac{3 EI}{L^3} & -\frac{3 EI}{L^2} &  0            & \frac{3 EI}{L^3}
\end{bmatrix}
\begin{bmatrix}
u_i\\
v_i\\
\theta_i\\
u_j\\
v_j\\
\end{bmatrix}
--->
![](https://i.upmath.me/svg/%5Crenewcommand%5Carraystretch%7B1.4%7D%0A%5Cbegin%7Bbmatrix%7D%0AX_i%5C%5C%0AY_i%5C%5C%0AM_i%5C%5C%0AX_j%5C%5C%0AY_j%5C%5C%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0A%5Cfrac%7BEA%7D%7BL%7D%20%26%20%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7BEA%7D%7BL%7D%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5C%5C%0A%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20%20%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%5C%5C%0A%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20%20%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B3%20EI%7D%7BL%7D%20%20%20%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%5C%5C%0A-%5Cfrac%7BEA%7D%7BL%7D%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7BEA%7D%7BL%7D%20%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5C%5C%0A%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%26%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Au_i%5C%5C%0Av_i%5C%5C%0A%5Ctheta_i%5C%5C%0Au_j%5C%5C%0Av_j%5C%5C%0A%5Cend%7Bbmatrix%7D)

y la matriz de rigidez para:

![](figs/Ke_rotula_empotrado.svg)

es:
<!---
Compile en: https://tex.s2cms.com

\renewcommand\arraystretch{1.4}
\begin{bmatrix}
X_i\\
Y_i\\
X_j\\
Y_j\\
M_j
\end{bmatrix}
=
\begin{bmatrix}
 \frac{EA}{L}  & 0                 & -\frac{EA}{L} & 0                 & 0                \\ 
  0            & \frac{3 EI}{L^3}  & 0             & -\frac{3 EI}{L^3} & \frac{3 EI}{L^2} \\ 
 -\frac{EA}{L} & 0                 & \frac{EA}{L}  & 0                 & 0                \\
  0            & -\frac{3 EI}{L^3} & 0             & \frac{3 EI}{L^3}  & -\frac{3 EI}{L^2}\\ 
  0            & \frac{3 EI}{L^2}  & 0             & -\frac{3 EI}{L^2} & \frac{3 EI}{L}  
\end{bmatrix}
\begin{bmatrix}
u_i\\
v_i\\
u_j\\
v_j\\
\theta_j
\end{bmatrix}
--->
![](https:////i.upmath.me/svg/%5Crenewcommand%5Carraystretch%7B1.4%7D%0A%5Cbegin%7Bbmatrix%7D%0AX_i%5C%5C%0AY_i%5C%5C%0AX_j%5C%5C%0AY_j%5C%5C%0AM_j%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0A%20%5Cfrac%7BEA%7D%7BL%7D%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5C%5C%20%0A%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%5C%5C%20%0A%20-%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7BEA%7D%7BL%7D%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5C%5C%0A%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E3%7D%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%5C%5C%20%0A%20%200%20%20%20%20%20%20%20%20%20%20%20%20%26%20%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%20%26%200%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20-%5Cfrac%7B3%20EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B3%20EI%7D%7BL%7D%20%20%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Au_i%5C%5C%0Av_i%5C%5C%0Au_j%5C%5C%0Av_j%5C%5C%0A%5Ctheta_j%0A%5Cend%7Bbmatrix%7D)

De otro lado, los programas:
* [ejemplo_empotrado_viga_rotula1.m](ejemplo_empotrado_viga_rotula1.m)
* [ejemplo_empotrado_viga_rotula2.m](ejemplo_empotrado_viga_rotula2.m)

ilustran la solución de vigas que tienen rótulas intermedias. 

**FALTA HACER LOS GRAFICOS DE LOS EJEMPLOS**