## Deducción de la matriz de rigidez de un pórtico en 2D

Los programas siguientes sirven para deducir la matriz de rigidez de un elemento finito de pórtico bidimensional:
* MATLAB: [deduccion_K_portico2D.m](deduccion_K_portico2D.m) (usa el toolbox de álgebra simbólica)
* PYTHON: [deduccion_K_portico2D.py](deduccion_K_portico2D.py)

<!---
Compile en: https://tex.s2cms.com

\renewcommand\arraystretch{1.4}
\begin{bmatrix}
X_i\\
Y_i\\
M_i\\
X_j\\
Y_j\\
M_j
\end{bmatrix}
=
\begin{bmatrix}
  \frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\
  0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\
  0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\
  -\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\
  0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\
  0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L}
\end{bmatrix}
\begin{bmatrix}
u_i\\
v_i\\
\theta_i\\
u_j\\
v_j\\
\theta_j
\end{bmatrix}
--->

![](https://tex.s2cms.ru/svg/%5Crenewcommand%5Carraystretch%7B1.4%7D%0A%5Cbegin%7Bbmatrix%7D%0AX_i%5C%5C%0AY_i%5C%5C%0AM_i%5C%5C%0AX_j%5C%5C%0AY_j%5C%5C%0AM_j%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0A%20%20%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%26%200%20%26%20-%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%26%200%20%5C%5C%0A%20%200%20%26%20%5Cfrac%7B12EI%7D%7BL%5E3%7D%20%26%20%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%200%20%26%20-%5Cfrac%7B12EI%7D%7BL%5E3%7D%20%26%20%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%5C%5C%0A%20%200%20%26%20%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B4EI%7D%7BL%7D%20%26%200%20%26%20-%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B2EI%7D%7BL%7D%20%5C%5C%0A%20%20-%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%26%200%20%26%20%5Cfrac%7BEA%7D%7BL%7D%20%26%200%20%26%200%20%5C%5C%0A%20%200%20%26%20-%5Cfrac%7B12EI%7D%7BL%5E3%7D%20%26%20-%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%200%20%26%20%5Cfrac%7B12EI%7D%7BL%5E3%7D%20%26%20-%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%5C%5C%0A%20%200%20%26%20%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B2EI%7D%7BL%7D%20%26%200%20%26%20-%5Cfrac%7B6EI%7D%7BL%5E2%7D%20%26%20%5Cfrac%7B4EI%7D%7BL%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Au_i%5C%5C%0Av_i%5C%5C%0A%5Ctheta_i%5C%5C%0Au_j%5C%5C%0Av_j%5C%5C%0A%5Ctheta_j%0A%5Cend%7Bbmatrix%7D)



## Cálculo de la carga nodal equivalente para una carga triangular
El programa [deduccion_fe_carga_triangular.m](deduccion_fe_carga_triangular.m) calcula el vector de fuerzas nodales equivalentes asociadas a la carga triangular mostrada en la figura:

![Image](fe_carga_triangular.svg)

para un elemento finito de pórtico bidimensional. Dicho programa utiliza el **toolbox de álgebra simbólica**.
