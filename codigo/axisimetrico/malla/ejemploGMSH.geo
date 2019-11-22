// se definen los puntos
B = 10;  // [m] ancho
H = 10;  // [m] alto
Point(1) = {0, -H/2, 0, 1.0};
Point(2) = {B, -H/2, 0, 1.0};
Point(3) = {B, +H/2, 0, 1.0};
Point(4) = {0, +H/2, 0, 0.1};

// se definen los bordes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// se define la superficie
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// se define la condición de borde
Physical Curve("empotrado") = {1, 2};
Physical Curve("eje") = {4};
Physical Curve("superior") = {3};

// se define la superficie
Physical Surface("suelo") = {1};

// Recombine the triangles into quads
Recombine Surface{1};

// Genere elementos serendipitos de 8 nodos
Mesh.SecondOrderIncomplete = 2;

// Removes all duplicate elementary entities (e.g., points having identical
// coordinates)
Coherence;

// Removes all duplicate mesh nodes.
Coherence Mesh;

// En GMSH hacer:
// 1. mesh->define->2D
// 2. mesh->define->Set order 2D
// 3. file->export->formato SU2