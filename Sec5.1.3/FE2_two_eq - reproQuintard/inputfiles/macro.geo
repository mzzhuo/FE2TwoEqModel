Mesh.MshFileVersion = 2.2;

L = 10; 
D = 0.3; 
lcar = 0.1; 

Point(1) = {0, 0, 0, lcar};
Point(2) = {L, 0, 0, lcar};
Point(3) = {L, D, 0, lcar};
Point(4) = {0, D, 0, lcar};

Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};


Line Loop(5) = {1, 4, -2, -3};
Plane Surface(6) = {5};

Transfinite Line { -1, -2 } = 41 Using Progression 0.98;
//Transfinite Line { -1, -2 } = 11 Using Bump .1; // 0.5 is better

Transfinite Line { 3,4 } = 2;

Transfinite Surface {6};
Recombine Surface {6};

Physical Surface("spe") = {6};


Physical Line("left")  = {3};
Physical Line("right") = {4};


Mesh 2;
//RefineMesh;