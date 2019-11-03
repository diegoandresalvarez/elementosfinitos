# Programa para extrapolar los esfuerzos en los EFs rectangulares de 8 y 9 nodos a partir de las estimaciones en los puntos de Gauss-Legendre

Para los esfuerzos calculados en los puntos I, II, III y IV y que se muestran a continuación:

![extrapolacion_esfuerzos.svg](extrapolacion_esfuerzos.svg | width=500)

la estimación de los esfuerzos en los nodos está dada por la siguiente fórmula de interpolación para el EF serendípito rectangular de 8 nodos:

La ecuación de interpolación buscada es:
<!---
Compile en: https://tex.s2cms.com

\begin{pmatrix}
\sigma_{\text{nodo }1} \\
\sigma_{\text{nodo }2} \\
\sigma_{\text{nodo }3} \\
\sigma_{\text{nodo }4} \\
\sigma_{\text{nodo }5} \\
\sigma_{\text{nodo }6} \\
\sigma_{\text{nodo }7} \\
\sigma_{\text{nodo }8}
\end{pmatrix}
 =
\underbrace{\begin{pmatrix} \frac{\sqrt{3}}{2} + 1 & - \frac{1}{2} & - \frac{1}{2} & 1 - \frac{\sqrt{3}}{2}\\ \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4}\\ - \frac{1}{2} & 1 - \frac{\sqrt{3}}{2} & \frac{\sqrt{3}}{2} + 1 & - \frac{1}{2}\\ \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4}\\ 1 - \frac{\sqrt{3}}{2} & - \frac{1}{2} & - \frac{1}{2} & \frac{\sqrt{3}}{2} + 1\\ \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4}\\ - \frac{1}{2} & \frac{\sqrt{3}}{2} + 1 & 1 - \frac{\sqrt{3}}{2} & - \frac{1}{2}\\ \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{\sqrt{3}}{4} + \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4} & \frac{1}{4} - \frac{\sqrt{3}}{4}
\end{pmatrix}}_{\ma{A}_2\ma{A}_1^{-1}}
\begin{pmatrix}
\sigma_I \\
\sigma_{II} \\
\sigma_{III} \\
\sigma_{IV} \\
\end{pmatrix}
--->
![](https://tex.s2cms.ru/svg/%5Cbegin%7Bpmatrix%7D%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D1%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D2%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D3%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D4%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D5%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D6%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D7%7D%20%5C%5C%0A%5Csigma_%7B%5Ctext%7Bnodo%20%7D8%7D%0A%5Cend%7Bpmatrix%7D%0A%20%3D%0A%5Cunderbrace%7B%5Cbegin%7Bpmatrix%7D%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%2B%201%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%201%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%5C%5C%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%5C%5C%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%201%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%2B%201%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%5C%5C%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%5C%5C%201%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%2B%201%5C%5C%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%5C%5C%20-%20%5Cfrac%7B1%7D%7B2%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%2B%201%20%26%201%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B2%7D%20%26%20-%20%5Cfrac%7B1%7D%7B2%7D%5C%5C%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%2B%20%5Cfrac%7B1%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%20%26%20%5Cfrac%7B1%7D%7B4%7D%20-%20%5Cfrac%7B%5Csqrt%7B3%7D%7D%7B4%7D%0A%5Cend%7Bpmatrix%7D%7D_%7B%5Cma%7BA%7D_2%5Cma%7BA%7D_1%5E%7B-1%7D%7D%0A%5Cbegin%7Bpmatrix%7D%0A%5Csigma_I%20%5C%5C%0A%5Csigma_%7BII%7D%20%5C%5C%0A%5Csigma_%7BIII%7D%20%5C%5C%0A%5Csigma_%7BIV%7D%20%5C%5C%0A%5Cend%7Bpmatrix%7D%5Cnonumber)

**Nota:** el orden I, II, III, IV corresponde al orden en que se evaluan los esfuerzos en el doble ciclo:
```matlab
i = 0;
for p = 1:n_gl
    for q = 1:n_gl
        i = i+1;
        % Se evalua el esfuerzo en el nodo i-ésimo de modo que:
        % i == 1 es equivalente a i = I
        % i == 2 es equivalente a i = II
        % i == 3 es equivalente a i = III
        % i == 4 es equivalente a i = IV              
    end
end
```

Por esta razón, la ecuación aquí deducida es diferente a la mostrada en Oñate (2009), página 329 a 334 y a la mostrada [aquí](http://books.google.com/books?id=lcSwbhop_XYC&pg=PA485&lpg=PA485&dq=%22nodal+stresses%22+%22gauss+points%22&source=bl&ots=75zUqMQDY1&sig=FJ_I-NbkkDkeKeIum9JOvlXqje4&hl=de&ei=M6OoTe_JJ-aJ0QH_haj5CA&sa=X&oi=book_result&ct=result&resnum=54&ved=0CJ4EEOgBMDU#v=onepage&q=%22nodal%20stresses%22%20%22gauss%20points%22&f=false).

Cödigo:
* MATLAB: [extrapolacion_esfuerzos_Q8_Q9.m](extrapolacion_esfuerzos_Q8_Q9.m)
* PYTHON: [extrapolacion_esfuerzos_Q8_Q9.py](extrapolacion_esfuerzos_Q8_Q9.py)