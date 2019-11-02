# Programa para extrapolar los esfuerzos en los EFs rectangulares de 8 y 9 nodos a partir de las estimaciones en los puntos de Gauss-Legendre

Para los esfuerzos calculados en los puntos I, II, III y IV y que se muestran a continuación:

![]()

la estimación de los esfuerzos en los nodos está dada por la siguiente fórmula de interpolación para el EF serendípito rectangular de 8 nodos:

La ecuación de interpolación buscada es:
<!---
Compile en: https://tex.s2cms.com

x = \frac{a+b}{2} + \frac{b-a}{2}\xi
--->
![](https://tex.s2cms.ru/svg/x%20%3D%20%5Cfrac%7Ba%2Bb%7D%7B2%7D%20%2B%20%5Cfrac%7Bb-a%7D%7B2%7D%5Cxi)

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