
//L = 1; 
L = 0.00357143; 
D = 0.00357143; 
//D = 0.0125; 
//D = 1; 

n_x = 1;
//n_x = 80;
//n_x = 280;          // no. of rve in x 
n_y = 1;          // no. of rve in y
n   = n_x * n_y;  // total no. of RVEs

lc1 = L/n_x/2/10;
lc2 = L/n_x/2/10;

Point(1) = {0.0, 0.0, 0.0, lc1};
Point(2) = {L/n_x  , 0.0, 0.0, lc1};
Point(3) = {L/n_x  , D/n_y  , 0.0, lc1};
Point(4) = {0.0, D/n_y  , 0.0, lc1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {1, 4};

Periodic Line {1,2} = {3,4};

linelooprve = newreg;
Line Loop(linelooprve) = {1,2,-3,-4};

x = L/n_x/2;
y = D/n_y/2;
r = L/n_x * 0.3250;

p1 = newp; Point(p1) = {x,   y,   0, lc2};
p2 = newp; Point(p2) = {x+r, y,   0, lc2};
p3 = newp; Point(p3) = {x,   y+r, 0, lc2};
p4 = newp; Point(p4) = {x-r, y,   0, lc2};
p5 = newp; Point(p5) = {x,   y-r, 0, lc2};

c1 = newl; Circle(c1) = {p2,p1,p3}; 
c2 = newl; Circle(c2) = {p3,p1,p4}; 
c3 = newl; Circle(c3) = {p4,p1,p5}; 
c4 = newl; Circle(c4) = {p5,p1,p2}; 

Periodic Line {c1} = {-c2};
Periodic Line {c1} = {-c4};
Periodic Line {c3} = {-c2};
Periodic Line {c3} = {-c4};

lineloopcir = newll;
Line Loop(lineloopcir) = {c1,c2,c3,c4};

psrve = newreg;
Plane Surface(psrve) = {linelooprve,lineloopcir};
Printf("Base Matrix Plane Surface: %g.", psrve);

pscircle = newreg;
Plane Surface(pscircle) = {lineloopcir};
Printf("Base Inclusion Plane Surface: %g.", pscircle);

Physical Surface("matrix") = {psrve};
Physical Surface("inclusion") = {pscircle};

//For i In {1:n_x-1}
//  my_new_surfs[] = Translate {i*L/n_x, 0, 0} { Duplicata { Surface{psrve,pscircle}; } };
//  Printf("‘%g'th: New surfaces ’%g’ and ’%g’", i, my_new_surfs[0], my_new_surfs[1]);
//  new_surfs_matrix[i-1] = my_new_surfs[0];
//  new_surfs_inclus[i-1] = my_new_surfs[1];
//EndFor

//Physical Surface("matrix") = {psrve,new_surfs_matrix[]};
//Physical Surface("inclusion") = {pscircle,new_surfs_inclus[]};

Physical Line("left")  = {4};
Physical Line("right") = {2};
//Physical Line("right") = {2795}; // manually pick

//Mesh 2;
//RefineMesh;


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// not copy rve and paste 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//L = 1; 
//D = 0.05; //0.00357143; 

//n_x = 20;          // no. of rve in x 
//n_y = 1;          // no. of rve in y
//n   = n_x * n_y;  // total no. of RVEs

//lc1 = L/n_x/5;
//lc2 = L/n_x/2/5;

//Point(1) = {0.0, 0.0, 0.0, lc1};
//Point(2) = {L  , 0.0, 0.0, lc1};
//Point(3) = {L  , D  , 0.0, lc1};
//Point(4) = {0.0, D  , 0.0, lc1};

//Line(1) = {1, 2};
//Line(2) = {4, 3};
//Line(3) = {1, 4};
//Line(4) = {2, 3};

//// define a Macro of circle for repeating
//Macro CheeseHole

//p1 = newp; Point(p1) = {x,   y,   0, lc2};
//p2 = newp; Point(p2) = {x+r, y,   0, lc2};
//p3 = newp; Point(p3) = {x,   y+r, 0, lc2};
//p4 = newp; Point(p4) = {x-r, y,   0, lc2};
//p5 = newp; Point(p5) = {x,   y-r, 0, lc2};

//c1 = newl; Circle(c1) = {p2,p1,p3}; 
//c2 = newl; Circle(c2) = {p3,p1,p4}; 
//c3 = newl; Circle(c3) = {p4,p1,p5}; 
//c4 = newl; Circle(c4) = {p5,p1,p2}; 

//lineloops[t] = newll;
//Line Loop(lineloops[t]) = {c1,c2,c3,c4};
//Printf("Line Loop number: lineloops[%g] = %g.", t, lineloops[t]);

//Return

//For j In {1:n_y}
//  For i In {1:n_x}
//  x = L/n_x/2 + (i - 1) * L/n_x;
//  y = D/n_y/2 + (j - 1) * D/n_y;
//  r = L/n_x * 0.3318;
//  t = i + (j - 1)*n_x - 1; // start from 0;
//  Call CheeseHole;
//  psloops[t] = newreg;
//  Plane Surface(psloops[t]) = {lineloops[t]};
//  EndFor
//EndFor

//rve = newll;
//Line Loop(rve) = {1, 4, -2, -3};

//Plane Surface(1) = {rve,lineloops[]};
//Physical Surface("matrix") = {1};

// //Plane Surface(2) = {lineloops[]};
// Physical Surface("inclusion") = {psloops[]};



//Physical Line("left")  = {3};
//Physical Line("right") = {4};
//Physical Line("lower") = {1};
//Physical Line("upper") = {2};

////Physical Point("corner") = {1,2,3,4};

////Mesh 2;
////RefineMesh;
