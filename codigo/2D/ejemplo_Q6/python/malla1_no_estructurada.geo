/*
  Generación de la malla de EFs para una viga utilizando EFs de 4 nodos 
  rectangulares
*/

// se definen los puntos
B = 0.8;   // [m] ancho
H = 0.28;  // [m] alto

// algunos parámetros de la malla
tm =  0.08;  // tamaño normal de la malla
tmr = 0.02;  // tamaño refinado de la malla

// se definen los puntos
Point(1) = {  0, 0, 0, tmr};
Point(2) = {  B, 0, 0, tmr};
Point(3) = {  B, H, 0, tm };
Point(4) = {0.6, H, 0, tmr};  // carga puntual
Point(5) = {0.4, H, 0, tmr};  // fin carga distribuída
Point(6) = {  0, H, 0, tm };

// se definen los bordes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// se define la superficie
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// se define la condición de borde
Physical Curve("carga_distribuida") = {5};

// se define la superficie
Physical Surface("solido") = {1};

// se definen los apoyos puntuales
Physical Point("apoyo_izq") = {1};
Physical Point("apoyo_der") = {2};

// se define el punto de la carga puntual
Physical Point("carga_puntual") = {4};

// use rectángulos en vez de triángulos
Recombine Surface{1};

// Se crea la malla
Mesh 2;

// Se usan rectángulos


// opciones de visualización
Mesh.SurfaceFaces = 1;
Mesh.Points = 1;

// se graba la malla (mesh output format)
Mesh.Format = 1;     Save "malla1.msh";    // gmsh
Mesh.Format = 39;    Save "malla1.inp";    // Abaqus
Mesh.Format = 42;    Save "malla1.su2";    // su2
Mesh.Format = 47;    Save "malla1.tochnog";
Mesh.Format = 50;    Save "malla1.matlab";
