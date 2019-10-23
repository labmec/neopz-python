// Gmsh project created on Mon Oct 21 14:10:39 2019
SetFactory("OpenCASCADE");

h=2.0;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {100.0, 0.0, 0.0, h};
Point(3) = {100.0, 10.0, 0.0, h};
Point(4) = {0.0, 10.0, 0.0, h};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

//+
Physical Surface("d_rock") = {1};
//+
Physical Curve("bc_nonFlux") = {4};
//+
Physical Curve("bc_nonFlux") += {2};
//+
Physical Curve("bc_Inlet") = {1};
//+
Physical Curve("bc_outlet") = {3};

//Coherence;

Transfinite Line {1,3} = 2 Using Bump 1;
Transfinite Surface "*";
Recombine Surface "*";
