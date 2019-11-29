// se definen los puntos
B = 0.8;   // [m] ancho
H = 0.28;  // [m] alto
Point(1) = {0, 0, 0, 0.01};
Point(2) = {B, 0, 0, 0.01};
Point(3) = {B, H, 0, 0.07};
Point(4) = {0.6, H, 0, 0.01};  // carga puntual
Point(5) = {0.4, H, 0, 0.01};  // fin carga distribuída
Point(6) = {0, H, 0, 0.07};

// se definen los bordes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// se define la superficie
Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// se define la condición de borde
Physical Curve("carga_distribuida") = {5};

// se define la superficie
Physical Surface("solido") = {1};

// Removes all duplicate elementary entities (e.g., points having identical
// coordinates)
Coherence;

// Removes all duplicate mesh nodes.
Coherence Mesh;

// En GMSH hacer:
// 1. mesh->define->2D
// 2. mesh->define->Set order 1D
// 3. file->export->formato SU2