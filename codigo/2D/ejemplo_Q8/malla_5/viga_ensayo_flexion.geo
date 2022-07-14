// Ejm 1: Malla no estructurada

tm = 0.03;
tmr = 0.01;

Point( 1) = {        0,    0, 0, tm};
Point( 2) = {0.10-0.01,    0, 0, tmr};
Point( 3) = {0.10+0.01,    0, 0, tmr};
Point( 4) = {1.30-0.01,    0, 0, tmr};
Point( 5) = {1.30+0.01,    0, 0, tmr};
Point( 6) = {     1.40,    0, 0, tm};
Point( 7) = {     1.40, 0.15, 0, tm};
Point( 8) = {0.90+0.01, 0.15, 0, tmr};
Point( 9) = {0.90-0.01, 0.15, 0, tmr};
Point(10) = {0.50+0.01, 0.15, 0, tmr};
Point(11) = {0.50-0.01, 0.15, 0, tmr};
Point(12) = {        0, 0.15, 0, tm};

Line( 1) = { 1,  2};
Line( 2) = { 2,  3};
Line( 3) = { 3,  4};
Line( 4) = { 4,  5};
Line( 5) = { 5,  6};
Line( 6) = { 6,  7};
Line( 7) = { 7,  8};
Line( 8) = { 8,  9};
Line( 9) = { 9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12,  1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

Plane Surface(1) = {1};

Physical Line("apoyos") = {2, 4};
Physical Line("cargas") = {8, 10};
Physical Surface("viga") = {1};

Recombine Surface{1}; // Para usar elementos cuadriláteros en vez de triángulos

// Ajustes finales de la malla		   
Mesh 2;                         // Generar la malla 2D
Mesh.SecondOrderIncomplete = 1; // Usar elementos finitos serendípitos (incompletos)
SetOrder 2;                     // Se definen elementos finitos de orden 2
Mesh.SurfaceFaces = 1;          // Ver las "caras" de los elementos finitos 2D
Mesh.Points = 1;                // Ver los nodos de la malla

Save "viga_ensayo_flexion.msh";
Save "viga_ensayo_flexion.m";