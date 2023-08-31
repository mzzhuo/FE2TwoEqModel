
//L = 0.00357143; 
//D = 0.00357143; 

L = 1/140; 
D = 1/140; 

//lc1 = L/20;
//lc2 = L/20;
lc1 = L/10;
lc2 = L/10;


Point(1) = {0.0, 0.0, 0.0, lc1};
Point(2) = {L  , 0.0, 0.0, lc1};
Point(3) = {L  , D  , 0.0, lc1};
Point(4) = {0.0, D  , 0.0, lc1};

Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};

// define a Macro of circle for repeating
Macro CheeseHole

p1 = newp; Point(p1) = {x,   y,   0, lc2};
p2 = newp; Point(p2) = {x+r, y,   0, lc2};
p3 = newp; Point(p3) = {x,   y+r, 0, lc2};
p4 = newp; Point(p4) = {x-r, y,   0, lc2};
p5 = newp; Point(p5) = {x,   y-r, 0, lc2};

c1 = newl; Circle(c1) = {p2,p1,p3}; 
c2 = newl; Circle(c2) = {p3,p1,p4}; 
c3 = newl; Circle(c3) = {p4,p1,p5}; 
c4 = newl; Circle(c4) = {p5,p1,p2}; 

lineloops[t] = newll;
Line Loop(lineloops[t]) = {c1,c2,c3,c4};
Printf("Line Loop number: lineloops[%g] = %g.", t, lineloops[t]);

Return

x = L/2; y = D/2;
r = L * 0.3250;
t = 0;
Call CheeseHole;
 psloops[t] = newreg;
 Plane Surface(psloops[t]) = {lineloops[t]};


rve = newll;
Line Loop(rve) = {1, 4, -2, -3};

Plane Surface(1) = {rve,lineloops[]};


Physical Surface("matrix") = {1};

Physical Surface("inclusion") = {psloops[]};

Physical Line("left")  = {3};
Physical Line("right") = {4};
Physical Line("lower") = {1};
Physical Line("upper") = {2};

Physical Point("corner") = {1,2,3,4};

Mesh 2;
//RefineMesh;