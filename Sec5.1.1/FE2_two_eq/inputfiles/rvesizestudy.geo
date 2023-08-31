

Point(1) = {0.00000000e+00,0.00000000e+00,0.0,3.57143000e-04};
Point(2) = {3.57143000e-03,0.00000000e+00,0.0,3.57143000e-04};
Point(3) = {3.57143000e-03,3.57143000e-03,0.0,3.57143000e-04};
Point(4) = {0.00000000e+00,3.57143000e-03,0.0,3.57143000e-04};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {4,3};
Line(4) = {1,4};


coordx[     1] = 2.05321592e-03; coordy[     1] = 1.61835663e-03;


Macro CheeseHole
p1 = newp; Point(p1) = {x,   y,   0, 3.57143000e-04};
p2 = newp; Point(p2) = {x+r, y,   0, 3.57143000e-04};
p3 = newp; Point(p3) = {x,   y+r, 0, 3.57143000e-04};
p4 = newp; Point(p4) = {x-r, y,   0, 3.57143000e-04};
p5 = newp; Point(p5) = {x,   y-r, 0, 3.57143000e-04};
c1 = newreg; Circle(c1) = {p2,p1,p3};
c2 = newreg; Circle(c2) = {p3,p1,p4};
c3 = newreg; Circle(c3) = {p4,p1,p5};
c4 = newreg; Circle(c4) = {p5,p1,p2};
lineloops[nll] = newreg;
Line Loop(lineloops[nll]) = {c1,c2,c3,c4};
Return


r = 1.16071475e-03;
n=      1;
For i In {1:n}
x = coordx[i];
y = coordy[i];
nll = i - 1;
Call CheeseHole;
intline[4*nll+0] = c1;
intline[4*nll+1] = c2;
intline[4*nll+2] = c3;
intline[4*nll+3] = c4;
EndFor
square = newreg;
Line Loop(square) = {1,2,-3,-4};
mat = newreg;
Plane Surface(mat) = {square,lineloops[]};
Physical Surface("matrix") = {mat};
Physical Line("int_mat") = {intline[]};


Delete lineloops;
Delete intline;
For i In {1:n}
x = coordx[i];
y = coordy[i];
nll = i - 1;
Call CheeseHole;
psloops[nll] = newreg;
Plane Surface(psloops[nll]) = {lineloops[nll]};
intline[4*nll+0] = c1;
intline[4*nll+1] = c2;
intline[4*nll+2] = c3;
intline[4*nll+3] = c4;
EndFor
Physical Surface("inclusion") = {psloops[]};
Physical Line("int_inc") = {intline[]};
Physical Line("left")  = {4};
Physical Line("right") = {2};
Physical Line("lower") = {1};
Physical Line("upper") = {3};
Physical Point("corner") = {1,2,3,4};
Mesh 2;
// RefineMesh;
