Point(0) = {-1,-1,0};
Point(1) = {1,-1,0};
Point(2) = {1,1,0};
Point(3) = {-1,1,0};


Line(4) = {0,1};
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,0};

Line Loop(1) = {4,5,6,7};
Surface(8) = {1};
Transfinite Surface{8};
Physical Surface(1) = {8};
Mesh(1);
Mesh(2);
Recombine Surface(8);